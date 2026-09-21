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
    args = p.parse_args()
    if not 0 <= args.pct_mito <= 100 or not 0 <= args.min_cells_fraction < 1:
        p.error('MT must be in [0,100], gene fraction in [0,1)')
    if min(args.n_pcs, args.n_neighbors, args.n_top_genes) < 2 or args.max_dense_gb <= 0:
        p.error('PCA/neighbors/HVG limits must be >=2 and memory limit positive')
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
