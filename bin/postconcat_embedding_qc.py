#!/usr/bin/env python3
"""Raw-count post-concatenation QC, with visualization-only Scanpy normalization.

No normalized matrix, PCA, or neighbor graph is written into the output MuData.
"""
import argparse
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def vector(x):
    return np.asarray(x).ravel()


def slug(x):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(x))


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)


def qc_observations(mdata, batch_key):
    gene = mdata.mod['gene']
    obs = gene.obs.copy()
    obs['total_counts'] = vector(gene.X.sum(axis=1))
    obs['n_genes_by_counts'] = vector((gene.X > 0).sum(axis=1))
    # These percentages are calculated before feature filtering upstream.
    # Recalculating after gene filtering would change the denominator.
    for key in ('pct_counts_mt', 'percent_mito'):
        if key in obs:
            obs['pct_counts_mt'] = pd.to_numeric(obs[key], errors='coerce')
            break
    else:
        raise ValueError('Missing pre-feature-filter mitochondrial percentage')
    if batch_key not in obs:
        if batch_key in mdata.obs:
            obs[batch_key] = mdata.obs[batch_key].reindex(obs.index)
        elif batch_key in mdata.mod['guide'].obs:
            obs[batch_key] = mdata.mod['guide'].obs[batch_key].reindex(obs.index)
        else:
            obs[batch_key] = 'all'
    assignment = mdata.mod['guide'].layers['guide_assignment']
    obs['n_guides_assigned'] = vector((assignment > 0).sum(axis=1))
    return obs


def measurement_panels(before, after, batch_key, cutoff, outdir):
    rows = []
    for number, label in enumerate(sorted(before[batch_key].astype(str).unique())):
        pre = before.loc[before[batch_key].astype(str) == label]
        post = after.loc[after[batch_key].astype(str) == label]
        fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex='col', sharey='col')
        for row, data, title in ((0, pre, 'Before mitochondrial filter'), (1, post, 'After mitochondrial filter')):
            ax = axes[row, 0]
            points = ax.scatter(np.log1p(data.total_counts), np.log1p(data.n_genes_by_counts),
                                c=data.pct_counts_mt, vmin=0, vmax=100, s=3, rasterized=True)
            fig.colorbar(points, ax=ax, label='Mitochondrial counts (%)')
            ax.set(xlabel='log1p RNA UMIs', ylabel='log1p detected genes', title=f'{title}: {len(data):,} cells')
            axes[row, 1].hist(data.pct_counts_mt, bins=np.linspace(0, 100, 51), color='#0284c7')
            axes[row, 1].axvline(cutoff, color='#dc2626', linestyle='--')
            axes[row, 1].set(xlabel='Mitochondrial counts (%)', ylabel='Cells')
            if len(data):
                axes[row, 2].boxplot([np.log1p(data.total_counts), np.log1p(data.n_genes_by_counts)],
                                    tick_labels=['RNA UMIs', 'Detected genes'], showfliers=False)
            axes[row, 2].set(ylabel='log1p value')
        fig.suptitle(f'{label}: post-concatenation QC, MT ≤ {cutoff:g}% | {len(pre):,} → {len(post):,} cells')
        save(fig, outdir / f'measurement_{number}_{slug(label)}_qc.png')
        rows.append(dict(measurement_set=label, cells_before=len(pre), cells_after=len(post),
                         cells_removed=len(pre)-len(post), threshold=cutoff, filter='pct_counts_mt <= threshold'))
    pd.DataFrame(rows).to_csv(outdir / 'measurement_filter_flow.tsv', sep='\t', index=False)


def gene_symbols(adata):
    """Return one symbol per feature without changing the delivered feature index."""
    for key in ('symbol', 'gene_symbol', 'gene_name'):
        if key in adata.var:
            values = adata.var[key].fillna('').astype(str).to_numpy()
            return np.where(values != '', values, adata.var_names.astype(str))
    return adata.var_names.astype(str).to_numpy()


def score_cell_cycle(temp, args):
    result = dict(requested=args.cell_cycle, reference=args.reference, tapseq_mode=args.tapseq_mode,
                  status='skipped')
    if args.cell_cycle == 'off':
        result['reason'] = 'disabled by QC_EMBEDDING_cell_cycle'
        return result
    if args.tapseq_mode and args.cell_cycle != 'on':
        result['reason'] = 'automatic scoring disabled for targeted TAP-seq panels'
        return result
    if args.reference not in ('human', 'mouse'):
        result['reason'] = f'no packaged markers for reference={args.reference}'
        return result
    markers = [line.strip() for line in args.cell_cycle_genes.read_text().splitlines() if line.strip()]
    if len(markers) < 97:
        raise ValueError(f'Expected 97 Regev cell-cycle markers, found {len(markers)}')
    s_markers, g2m_markers = markers[:43], markers[43:]
    if args.reference == 'mouse':
        s_markers = [gene.title() for gene in s_markers]
        g2m_markers = [gene.title() for gene in g2m_markers]
        result['marker_conversion'] = 'human symbols converted to mouse-style capitalization'
    symbols = pd.Index(gene_symbols(temp))
    # A dedicated temporary view lets Scanpy match symbols while the pipeline's
    # raw-count MuData retains stable Ensembl feature identifiers.
    score_names = symbols.where(symbols.notna() & (symbols != ''), temp.var_names.astype(str))
    score_names = pd.Index(score_names).astype(str)
    score_view = ad.AnnData(X=temp.X, obs=temp.obs.copy(), var=pd.DataFrame(index=score_names))
    score_view.var_names_make_unique()
    s_present = [gene for gene in s_markers if gene in score_view.var_names]
    g2m_present = [gene for gene in g2m_markers if gene in score_view.var_names]
    result.update(marker_source=str(args.cell_cycle_genes), s_markers_total=len(s_markers),
                  g2m_markers_total=len(g2m_markers), s_markers_present=len(s_present),
                  g2m_markers_present=len(g2m_present))
    if min(len(s_present), len(g2m_present)) < args.min_cell_cycle_genes:
        result['reason'] = (f'only {len(s_present)} S and {len(g2m_present)} G2M markers present; '
                            f'minimum is {args.min_cell_cycle_genes} per phase')
        return result
    sc.tl.score_genes_cell_cycle(score_view, s_genes=s_present, g2m_genes=g2m_present,
                                 use_raw=False, random_state=args.seed)
    for key in ('S_score', 'G2M_score', 'phase'):
        temp.obs[key] = score_view.obs[key].to_numpy()
    counts = temp.obs['phase'].value_counts().reindex(['G1', 'S', 'G2M'], fill_value=0)
    pd.DataFrame({'phase': counts.index, 'cells': counts.values,
                  'fraction': counts.values / temp.n_obs}).to_csv(
                      args.outdir / 'cell_cycle_phase_counts.tsv', sep='\t', index=False)
    table = pd.crosstab(temp.obs[args.batch_key].astype(str), temp.obs['phase'], normalize='index')
    table = table.reindex(columns=['G1', 'S', 'G2M'], fill_value=0)
    ax = table.plot.bar(stacked=True, figsize=(max(8, .45 * len(table)), 5),
                        color=['#64748b', '#0ea5e9', '#f97316'])
    ax.set(xlabel='Measurement set', ylabel='Fraction of cells', ylim=(0, 1),
           title=f'{args.stage}: cell-cycle phase by measurement set')
    ax.legend(title='Phase', bbox_to_anchor=(1.01, 1), loc='upper left')
    save(ax.figure, args.outdir / 'cell_cycle_by_measurement_set.png')
    result.update(status='completed', phase_counts={str(k): int(v) for k, v in counts.items()})
    return result


def faceted_measurement_pca(temp, args):
    labels = sorted(temp.obs[args.batch_key].astype(str).unique())
    colors = plt.get_cmap('tab20')(np.linspace(0, 1, max(2, len(labels))))
    palette = dict(zip(labels, colors))
    cols = min(5, max(1, len(labels)))
    rows = int(np.ceil(len(labels) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.6 * rows), squeeze=False,
                             sharex=True, sharey=True)
    xy = temp.obsm['X_pca'][:, :2]
    values = temp.obs[args.batch_key].astype(str).to_numpy()
    for ax, label in zip(axes.ravel(), labels):
        selected = values == label
        ax.scatter(xy[:, 0], xy[:, 1], s=.35, c='#cbd5e1', alpha=.12, linewidths=0, rasterized=True)
        ax.scatter(xy[selected, 0], xy[selected, 1], s=.8, color=palette[label], alpha=.6,
                   linewidths=0, rasterized=True)
        ax.set_title(f'{label}\nn={selected.sum():,}')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
    for ax in axes.ravel()[len(labels):]:
        ax.axis('off')
    fig.suptitle(f'{args.stage}: global PCA highlighted by measurement set')
    save(fig, args.outdir / 'pca_by_measurement_set.png')
    # Dedicated pooled view in the visual style of the reference QC: every
    # retained high-quality cell is shown and colored by measurement set.
    fig, ax = plt.subplots(figsize=(12, 8))
    order = np.random.default_rng(args.seed).permutation(temp.n_obs)
    for label in labels:
        selected = order[values[order] == label]
        ax.scatter(xy[selected, 0], xy[selected, 1], s=1.2, color=palette[label], alpha=.48,
                   linewidths=0, rasterized=True, label=f'{label} (n={len(selected):,})')
    ax.set(xlabel='PC1', ylabel='PC2',
           title=f'{args.stage}: all {temp.n_obs:,} high-quality cells colored by measurement set')
    ax.legend(title='Measurement set', bbox_to_anchor=(1.01, 1), loc='upper left',
              markerscale=5, frameon=False, fontsize=8)
    save(fig, args.outdir / 'pca_measurement_sets_colored.png')


def run_leiden_sweep(temp, args):
    result = dict(enabled=args.enable_leiden, requested_resolutions=args.leiden_resolutions,
                  diagnostic_resolution=args.leiden_diagnostic_resolution,
                  n_iterations=args.leiden_n_iterations, implementation='python-igraph community_leiden')
    if not args.enable_leiden:
        result.update(status='skipped', reason='disabled')
        return result
    import igraph as ig
    adjacency = sp.csr_matrix(temp.obsp['connectivities'])
    graph = ig.Graph.Weighted_Adjacency(adjacency, mode='undirected', attr='weight', loops=False)
    rows = []
    keys = []
    for resolution in args.leiden_resolutions:
        key = f'leiden_res_{resolution:g}'
        partition = graph.community_leiden(objective_function='modularity', weights='weight',
                                           resolution=resolution, n_iterations=args.leiden_n_iterations,
                                           beta=0.01, initial_membership=None)
        labels = pd.Categorical(np.asarray(partition.membership).astype(str))
        temp.obs[key] = labels
        sizes = pd.Series(labels).value_counts()
        rows.append(dict(resolution=resolution, clusters=len(sizes), smallest_cluster=int(sizes.min()),
                         median_cluster=float(sizes.median()), largest_cluster=int(sizes.max()),
                         modularity=float(partition.modularity)))
        keys.append(key)
    pd.DataFrame(rows).to_csv(args.outdir / 'leiden_resolution_summary.tsv', sep='\t', index=False)
    temp.obs[[args.batch_key] + keys].to_csv(args.outdir / 'leiden_assignments.tsv.gz', sep='\t')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sweep = pd.DataFrame(rows)
    axes[0].plot(sweep.resolution, sweep.clusters, 'o-', color='#0284c7')
    axes[0].set(xlabel='Leiden resolution', ylabel='Number of clusters')
    axes[1].plot(sweep.resolution, sweep.modularity, 'o-', color='#7c3aed')
    axes[1].set(xlabel='Leiden resolution', ylabel='Modularity')
    save(fig, args.outdir / 'leiden_resolution_sweep.png')
    cols = min(4, len(keys))
    fig, axes = plt.subplots(int(np.ceil(len(keys)/cols)), cols,
                             figsize=(5*cols, 4*int(np.ceil(len(keys)/cols))), squeeze=False)
    for ax, key, resolution in zip(axes.ravel(), keys, args.leiden_resolutions):
        sc.pl.umap(temp, color=key, ax=ax, show=False, size=4, legend_loc='none',
                   title=f'Leiden resolution {resolution:g}')
    for ax in axes.ravel()[len(keys):]:
        ax.axis('off')
    fig.suptitle(f'{args.stage}: Leiden resolution sweep')
    save(fig, args.outdir / 'leiden_sweep_umap.png')
    diagnostic_index = int(np.argmin(np.abs(np.asarray(args.leiden_resolutions) - args.leiden_diagnostic_resolution)))
    diagnostic_key = keys[diagnostic_index]
    composition = pd.crosstab(temp.obs[diagnostic_key], temp.obs[args.batch_key].astype(str), normalize='index')
    ax = composition.plot.bar(stacked=True, figsize=(max(9, .5*len(composition)), 5), colormap='tab20')
    ax.set(xlabel=f'Leiden cluster ({diagnostic_key})', ylabel='Measurement-set fraction', ylim=(0, 1),
           title=f'{args.stage}: measurement-set composition by Leiden cluster')
    ax.legend(title=args.batch_key, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=7)
    save(ax.figure, args.outdir / 'leiden_measurement_set_composition.png')
    result.update(status='completed', effective_diagnostic_resolution=args.leiden_resolutions[diagnostic_index],
                  diagnostic_key=diagnostic_key, diagnostic_clusters=int(temp.obs[diagnostic_key].nunique()))
    return result


def embedding_plots(raw, obs, args):
    # Construct an independent temporary object, not a view of delivered counts.
    temp = ad.AnnData(X=sp.csr_matrix(raw.X, dtype=np.float32).copy(), obs=obs.copy(), var=raw.var.copy())
    temp = temp[vector(temp.X.sum(1)) > 0, vector((temp.X > 0).sum(0)) > 0].copy()
    result = dict(cells=temp.n_obs, genes=temp.n_vars, requested_pcs=args.n_pcs,
                  requested_neighbors=args.n_neighbors, seed=args.seed,
                  normalization='normalize_total(target_sum=None); log1p', batch_correction=False)
    if min(temp.shape) < 3:
        result.update(status='skipped', reason='Fewer than three nonzero cells or genes')
        return result
    result['normalization_target_sum'] = float(np.median(vector(temp.X.sum(1))))
    original_totals = vector(temp.X.sum(1))
    sc.pp.normalize_total(temp, target_sum=None)
    sc.pp.log1p(temp)
    result['cell_cycle'] = score_cell_cycle(temp, args)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(original_totals, bins=60)
    axes[0].set(title='Raw RNA counts', xlabel='UMIs per cell', ylabel='Cells')
    axes[1].hist(vector(temp.X.sum(1)), bins=60)
    axes[1].set(title='Median-depth normalization + log1p', xlabel='Sum of log-normalized expression', ylabel='Cells')
    save(fig, args.outdir / 'normalization_check.png')
    if args.hvg_batch_key and args.hvg_batch_key not in temp.obs:
        raise ValueError(f'HVG batch key is unavailable: {args.hvg_batch_key}')
    if temp.n_vars > args.n_top_genes:
        sc.pp.highly_variable_genes(temp, flavor='seurat', n_top_genes=args.n_top_genes,
                                   batch_key=args.hvg_batch_key or None)
        temp.var.to_csv(args.outdir / 'hvg_selection.tsv', sep='\t')
        temp = temp[:, temp.var.highly_variable].copy()
        result['feature_selection'] = 'seurat HVG'
    else:
        result['feature_selection'] = 'all expressed genes (panel smaller than HVG limit)'
    # Constant genes cannot contribute to PCA and can cause numerical failures.
    means = vector(temp.X.mean(0))
    variance = vector(temp.X.multiply(temp.X).mean(0)) - means ** 2
    temp = temp[:, variance > 1e-8].copy()
    pcs = min(args.n_pcs, temp.n_obs - 1, temp.n_vars - 1)
    if pcs < 2:
        result.update(status='skipped', reason='Fewer than two usable PCA dimensions')
        return result
    # Scaling densifies HVGs, so fail clearly instead of silently subsampling cells.
    required_gb = temp.n_obs * temp.n_vars * 4 / 1e9
    if required_gb > args.max_dense_gb:
        raise ValueError(f'Scaled HVG array requires {required_gb:.2f} GB; increase QC_EMBEDDING_max_dense_gb and task memory')
    sc.pp.scale(temp, max_value=10)
    sc.pp.pca(temp, n_comps=pcs, svd_solver='arpack', random_state=args.seed)
    faceted_measurement_pca(temp, args)
    vr = temp.uns['pca']['variance_ratio']
    pd.DataFrame({'PC': np.arange(1, len(vr)+1), 'variance_ratio': vr,
                  'cumulative_variance': np.cumsum(vr)}).to_csv(args.outdir / 'pca_variance_ratio.tsv', sep='\t', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].semilogy(np.arange(1, len(vr)+1), vr, 'o-')
    axes[0].set(xlabel='PC', ylabel='Variance explained (log scale)')
    axes[1].plot(np.arange(1, len(vr)+1), np.cumsum(vr), 'o-')
    axes[1].set(xlabel='PC', ylabel='Cumulative variance explained')
    save(fig, args.outdir / 'pca_variance_ratio.png')
    neighbors = min(args.n_neighbors, temp.n_obs-1)
    sc.pp.neighbors(temp, n_neighbors=neighbors, n_pcs=pcs, random_state=args.seed)
    sc.tl.umap(temp, random_state=args.seed)
    result['leiden'] = run_leiden_sweep(temp, args)
    colors = list(dict.fromkeys([args.batch_key, 'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'n_guides_assigned'] +
                 [k for k in ('pct_counts_ribo', 'doublet_score', 'rep', 'well', 'diff_day', 'guide_category', 'S_score', 'G2M_score') if k in temp.obs]))
    for basis in ('pca', 'umap'):
        fig, axes = plt.subplots(int(np.ceil(len(colors)/3)), 3, figsize=(18, 5*int(np.ceil(len(colors)/3))))
        for ax, key in zip(np.asarray(axes).ravel(), colors):
            if not pd.api.types.is_numeric_dtype(temp.obs[key]):
                temp.obs[key] = temp.obs[key].astype(str).astype('category')
            sc.pl.embedding(temp, basis=basis, color=key, ax=ax, show=False, size=5,
                            legend_loc='right margin' if temp.obs[key].nunique() <= 20 else 'none')
        for ax in np.asarray(axes).ravel()[len(colors):]:
            ax.axis('off')
        fig.suptitle(f'{args.stage}: {temp.n_obs:,} cells; {pcs} PCs; {neighbors} neighbors; no batch correction')
        save(fig, args.outdir / f'{basis}_qc_panel.png')
    result.update(status='completed', effective_pcs=pcs, effective_neighbors=neighbors, pca_features=temp.n_vars,
                  omitted_covariates=[k for k in ('pct_counts_ribo', 'S_score', 'G2M_score') if k not in temp.obs])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input_mudata')
    p.add_argument('--output', default='postconcat_filtered.h5mu')
    p.add_argument('--outdir', type=Path, default=Path('embedding_qc'))
    p.add_argument('--stage', default='before_clone')
    p.add_argument('--pct-mito', type=float, default=25)
    p.add_argument('--min-cells-fraction', type=float, default=0.05)
    p.add_argument('--n-pcs', type=int, default=50)
    p.add_argument('--n-neighbors', type=int, default=15)
    p.add_argument('--n-top-genes', type=int, default=6000)
    p.add_argument('--hvg-batch-key', default='')
    p.add_argument('--batch-key', default='batch')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--max-dense-gb', type=float, default=16)
    p.add_argument('--reference', choices=('human', 'mouse'), default='human')
    p.add_argument('--tapseq-mode', action=argparse.BooleanOptionalAction, default=False)
    p.add_argument('--cell-cycle', choices=('auto', 'on', 'off'), default='auto')
    p.add_argument('--cell-cycle-genes', type=Path, required=True)
    p.add_argument('--min-cell-cycle-genes', type=int, default=10)
    p.add_argument('--enable-leiden', action=argparse.BooleanOptionalAction, default=True)
    p.add_argument('--leiden-resolutions', default='0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2')
    p.add_argument('--leiden-diagnostic-resolution', type=float, default=0.5)
    p.add_argument('--leiden-n-iterations', type=int, default=2)
    args = p.parse_args()
    try:
        args.leiden_resolutions = [float(value) for value in args.leiden_resolutions.split(',') if value.strip()]
    except ValueError as error:
        p.error(f'Invalid Leiden resolution list: {error}')
    if not 0 <= args.pct_mito <= 100 or not 0 <= args.min_cells_fraction < 1:
        p.error('MT must be in [0,100], gene fraction in [0,1)')
    if min(args.n_pcs, args.n_neighbors, args.n_top_genes) < 2 or args.max_dense_gb <= 0:
        p.error('PCA/neighbors/HVG limits must be >=2 and memory limit positive')
    if args.min_cell_cycle_genes < 1 or args.leiden_n_iterations < 1 or not args.leiden_resolutions:
        p.error('Cell-cycle minimum and Leiden iterations/resolution list must be positive/nonempty')
    if min(args.leiden_resolutions) < 0 or args.leiden_diagnostic_resolution < 0:
        p.error('Leiden resolutions must be nonnegative')
    args.outdir.mkdir(parents=True, exist_ok=True)
    mdata = md.read_h5mu(args.input_mudata)
    for mod in mdata.mod.values():
        if not mod.obs_names.equals(mdata.obs_names):
            raise ValueError('Modalities must share the same ordered cell intersection')
    obs = qc_observations(mdata, args.batch_key)
    keep = np.isfinite(obs.pct_counts_mt) & (obs.pct_counts_mt <= args.pct_mito)
    measurement_panels(obs, obs.loc[keep], args.batch_key, args.pct_mito, args.outdir)
    filtered = mdata[keep.to_numpy()].copy()
    if filtered.n_obs == 0:
        raise ValueError('No cells remain after mitochondrial QC')
    summary = embedding_plots(filtered.mod['gene'], obs.loc[keep], args)
    gene = filtered.mod['gene']
    support = vector((gene.X > 0).sum(0))
    required = max(1, int(np.floor(filtered.n_obs * args.min_cells_fraction)) + 1)
    summary.update(stage=args.stage, input_cells=mdata.n_obs, retained_cells=filtered.n_obs,
                   mito_cutoff=args.pct_mito, genes_before=gene.n_vars, minimum_gene_cells=required,
                   genes_after=int((support >= required).sum()), normalized_matrix_saved=False)
    filtered.mod['gene'] = gene[:, support >= required].copy()
    if filtered.mod['gene'].n_vars == 0:
        raise ValueError('No genes remain after fractional gene-support QC')
    from concat_preprocessed_rna import recompute_gene_metrics
    recompute_gene_metrics(filtered.mod['gene'])
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis('off')
    flow = [f'{args.stage}: {mdata.n_obs:,} qualified raw-count cells',
            f'MT ≤ {args.pct_mito:g}% → {filtered.n_obs:,} cells',
            f'Temporary median-depth normalization → log1p → PCA → UMAP\n{summary["status"]}; no normalized matrix saved',
            f'Raw counts: gene support ≥ {required:,} cells\n{gene.n_vars:,} → {filtered.mod["gene"].n_vars:,} genes']
    for index, label in enumerate(flow):
        y = 0.9 - index * 0.24
        ax.text(.5, y, label, ha='center', va='center', fontsize=11,
                bbox=dict(boxstyle='round,pad=.6', facecolor='#f0f9ff', edgecolor='#0284c7'))
        if index:
            ax.annotate('', xy=(.5, y+.075), xytext=(.5, y+.16), arrowprops=dict(arrowstyle='->'))
    save(fig, args.outdir / 'postconcat_qc_flow.png')
    filtered.update()
    filtered.write_h5mu(args.output)
    (args.outdir / 'embedding_qc_metrics.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
