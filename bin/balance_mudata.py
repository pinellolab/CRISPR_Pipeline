#!/usr/bin/env python3
"""Exact cell and full-transcriptome RNA-UMI matching; preserve guide calls."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from qc_mudata_io import read_mudata_without_uns


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def axis(group):
    return pd.Index(group[group.attrs['_index']].asstr()[:])


def read_csr(group):
    return sparse.csr_matrix((group['data'][:], group['indices'][:], group['indptr'][:]),
                             shape=group.attrs['shape'])


def select_cells(cells, number, seed):
    cells = pd.Index(cells)
    if not cells.is_unique or not 0 < number <= len(cells):
        raise ValueError('Cell IDs must be unique and requested count attainable')
    return sorted(cells, key=lambda c: (hashlib.sha256(f'{seed}\0{c}'.encode()).digest(), c))[:number]


def canonical_counts(matrix):
    result = sparse.csr_matrix(matrix, copy=True)
    result.sum_duplicates(); result.sort_indices(); result.eliminate_zeros()
    values = result.data
    if (not np.isfinite(values).all() or np.any(values < 0)
            or not np.array_equal(values, np.rint(values))):
        raise ValueError('RNA X must contain nonnegative finite integer raw counts')
    return result


def exact_thin(matrix, budget, seed):
    matrix = canonical_counts(matrix)
    original = matrix.data.astype(np.int64)
    total = int(original.sum(dtype=np.int64))
    if not 0 < budget <= total:
        raise ValueError(f'RNA budget {budget} is not attainable from {total}')
    if budget != total:
        draw = np.random.Generator(np.random.PCG64(seed)).multivariate_hypergeometric(
            original, int(budget), method='count')
        if int(draw.sum(dtype=np.int64)) != budget or np.any(draw < 0) or np.any(draw > original):
            raise RuntimeError('Hypergeometric sampling invariants failed')
        matrix.data = draw.astype(np.int32 if draw.max(initial=0) <= np.iinfo(np.int32).max else np.int64)
        matrix.eliminate_zeros()
    if int(matrix.sum(dtype=np.int64)) != budget:
        raise RuntimeError('Sampled matrix RNA total differs from target')
    return matrix, total


def signature(matrix):
    matrix = canonical_counts(matrix)
    digest = hashlib.sha256()
    for values in (np.asarray(matrix.shape, dtype=np.int64), matrix.indptr.astype(np.int64),
                   matrix.indices.astype(np.int64), matrix.data.astype(np.int64)):
        digest.update(values.tobytes())
    return digest.hexdigest()


def update_rna(gene):
    x = sparse.csr_matrix(gene.X)
    total = np.asarray(x.sum(1)).ravel().astype(np.int64)
    detected = np.diff(x.indptr)
    for key in ('total_gene_umis', 'total_counts', 'n_counts'):
        gene.obs[key] = total
    for key in ('num_expressed_genes', 'n_genes_by_counts'):
        gene.obs[key] = detected
    for key in ('log1p_total_counts', 'log_total_gene_umis'):
        gene.obs[key] = np.log1p(total)
    for key in ('log1p_n_genes_by_counts', 'log_num_expressed_genes'):
        gene.obs[key] = np.log1p(detected)
    column = next((key for key in ('gene_name', 'gene_symbol', 'feature_name') if key in gene.var), None)
    symbols = gene.var[column].fillna('').astype(str) if column else pd.Series(gene.var_names, index=gene.var_names)
    for label, mask in (('mt', symbols.str.upper().str.startswith('MT-').to_numpy()),
                        ('ribo', symbols.str.upper().str.match(r'^RP[SL]').to_numpy())):
        counts = np.asarray(x[:, mask].sum(1)).ravel()
        pct = np.divide(100 * counts, total, out=np.zeros(len(total), dtype=float), where=total > 0)
        gene.obs[f'total_counts_{label}'] = counts
        gene.obs[f'log1p_total_counts_{label}'] = np.log1p(counts)
        gene.obs[f'pct_counts_{label}'] = pct
        if label == 'mt': gene.obs['percent_mito'] = pct
    sums = np.asarray(x.sum(0)).ravel()
    gene.var['total_counts'] = sums
    gene.var['n_cells_by_counts'] = np.asarray((x > 0).sum(0)).ravel()
    gene.var['mean_counts'] = sums / gene.n_obs
    gene.var['log1p_mean_counts'] = np.log1p(sums / gene.n_obs)
    gene.var['log1p_total_counts'] = np.log1p(sums)
    gene.var['pct_dropout_by_counts'] = 100 * (1 - gene.var['n_cells_by_counts'] / gene.n_obs)
    # Recent AnnData versions expose X as the None layer key. Preserve it.
    for key in list(gene.layers):
        if key is not None:
            del gene.layers[key]
    gene.obsm.clear(); gene.obsp.clear(); gene.varm.clear(); gene.varp.clear()
    gene.uns.clear(); gene.raw = None
    if gene.X is None:
        raise RuntimeError('Raw RNA matrix was lost during derived-state cleanup')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True, type=Path)
    p.add_argument('--parent-inference', required=True, type=Path)
    p.add_argument('--dataset', required=True)
    p.add_argument('--cells', required=True, type=int)
    p.add_argument('--rna-umis', required=True, type=int)
    p.add_argument('--seed', required=True, type=int)
    p.add_argument('--outdir', type=Path, default=Path('.'))
    args = p.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    outputs = [args.outdir / name for name in ('sampled_full_mudata.h5mu', 'balanced_inference_input.h5mu', 'downsampling_report.json')]
    if any(path.exists() for path in outputs):
        raise ValueError('Refusing to overwrite existing balanced outputs; use a new directory')
    with h5py.File(args.parent_inference, 'r') as handle:
        eligible = axis(handle['mod/gene/obs'])
        tested_genes = axis(handle['mod/gene/var'])
        parent_guides = axis(handle['mod/guide/var'])
    source = read_mudata_without_uns(args.source)
    selected = select_cells(eligible, args.cells, args.seed)
    for modality in source.mod.values():
        if not modality.obs_names.is_unique or not pd.Index(selected).isin(modality.obs_names).all():
            raise ValueError('Selected cell IDs are not present uniquely in every modality')
    sample = source[selected, :].copy()
    if not sample['guide'].var_names.equals(parent_guides):
        raise ValueError('Parent inference and source guide feature identities differ')
    if 'guide_assignment' not in sample['guide'].layers:
        raise ValueError('Source lacks guide_assignment')
    guide_counts = signature(sample['guide'].X)
    guide_calls = signature(sample['guide'].layers['guide_assignment'])
    with h5py.File(args.parent_inference, 'r') as handle:
        guide = handle['mod/guide']
        idx = axis(guide['obs']).get_indexer(selected)
        if np.any(idx < 0):
            raise ValueError('Selected cells missing from parent guide modality')
        if (guide_counts != signature(read_csr(guide['X'])[idx])
                or guide_calls != signature(read_csr(guide['layers/guide_assignment'])[idx])):
            raise ValueError('Parent and source guide counts or assignments differ')
    sample['gene'].X, before = exact_thin(sample['gene'].X, args.rna_umis, args.seed)
    update_rna(sample['gene'])
    sample.uns.clear()
    for key in list(sample.obs):
        if key.startswith('gene:') and key[5:] in sample['gene'].obs:
            sample.obs[key] = sample['gene'].obs[key[5:]].to_numpy()
    sample.update()
    if sample.n_obs != args.cells or not sample['gene'].obs_names.equals(pd.Index(selected)):
        raise RuntimeError('MuData cell order/count changed during sampling')
    mask = sample['gene'].var_names.isin(tested_genes)
    if int(mask.sum()) != len(tested_genes):
        raise ValueError('Parent analysis gene universe is not contained in full source')
    analysis = sample.copy()
    analysis.mod['gene'] = sample['gene'][:, mask].copy()
    analysis.update()
    if guide_counts != signature(analysis['guide'].X) or guide_calls != signature(analysis['guide'].layers['guide_assignment']):
        raise RuntimeError('Guide counts or assignments changed')
    report = dict(dataset=args.dataset, seed=args.seed, numpy_version=np.__version__,
                  rna_sampling_method='multivariate_hypergeometric', algorithm='count', bit_generator='PCG64',
                  source=str(args.source.resolve()), parent_inference=str(args.parent_inference.resolve()),
                  source_sha256=sha256(args.source), parent_inference_sha256=sha256(args.parent_inference),
                  eligible_cells=len(eligible), selected_cells=args.cells,
                  selected_cells_sha256=hashlib.sha256('\n'.join(selected).encode()).hexdigest(),
                  source_selected_rna_umis=before, target_total_rna_umis=args.rna_umis,
                  observed_total_rna_umis=int(sample['gene'].X.sum(dtype=np.int64)),
                  mean_full_transcriptome_rna_umis=args.rna_umis / args.cells,
                  full_transcriptome_genes=sample['gene'].n_vars, analysis_genes=analysis['gene'].n_vars,
                  analysis_gene_policy='freeze parent full-depth tested gene universe; no new cell filtering',
                  analysis_gene_ids_sha256=hashlib.sha256('\n'.join(analysis['gene'].var_names).encode()).hexdigest(),
                  analysis_gene_rna_umis=int(analysis['gene'].X.sum(dtype=np.int64)),
                  zero_rna_cells=int((sample['gene'].obs.total_gene_umis == 0).sum()),
                  guide_counts_sha256=guide_counts, guide_assignment_sha256=guide_calls,
                  guide_counts_unchanged=True, guide_assignment_unchanged=True,
                  full_rna_matrix_sha256=signature(sample['gene'].X),
                  source_qc='Inherited parent cell/guide eligibility; no remapping, doublet calling or guide reassignment',
                  comparison_caveat='Equal cell count and global RNA depth do not equalize biology, per-cell depth distributions or assay-specific gene detection.')
    sample.uns['balanced_matching'] = report.copy(); analysis.uns['balanced_matching'] = report.copy()
    sample.write_h5mu(outputs[0]); analysis.write_h5mu(outputs[1])
    pd.DataFrame({'cell_id': selected}).to_csv(args.outdir / 'selected_cells.tsv', sep='\t', index=False)
    outputs[2].write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__': main()
