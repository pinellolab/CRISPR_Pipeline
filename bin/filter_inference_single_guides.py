#!/usr/bin/env python3
"""Select single-guide cells for low-MOI inference after high-MOI assignment, without refitting calls or modifying counts.

A cell is kept when (1) the binary guide_assignment layer calls exactly one guide, (2) that guide carries at least --min-umis
raw guide UMIs, and (3) it carries at least --min-ratio times the UMIs of the largest other guide in the cell. The ratio is what
separates a genuine second infection from the diffuse one-UMI background seen in low-depth guide libraries: a second guide at a
quarter of a 600-UMI top guide fails it, fourteen one-UMI background guides do not. With --min-umis 0 --min-ratio 0 the selection
reduces to "exactly one assigned guide" (the earlier behaviour). Raw counts and assignment layers are never changed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
try:
    from anndata.io import read_elem
except ImportError:  # AnnData 0.10 in some pinned pipeline environments.
    from anndata._io.specs import read_elem
from scipy import sparse

from qc_mudata_io import read_mudata_without_uns


def signature(matrix):
    """SHA-256 of a count matrix in canonical CSR form, so a report can prove the counts left the step untouched."""
    matrix = sparse.csr_matrix(matrix, copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    digest = hashlib.sha256()
    for values in (np.asarray(matrix.shape, dtype=np.int64), matrix.indptr.astype(np.int64),
                   matrix.indices.astype(np.int64), matrix.data.astype(np.float64)):
        digest.update(values.tobytes())
    return digest.hexdigest()


def select_single_guide_cells(data, min_umis=5, min_ratio=3.0):
    if not {'gene', 'guide'}.issubset(data.mod):
        raise ValueError('gene and guide modalities are required')
    if min_umis < 0 or min_ratio < 0:
        raise ValueError('min_umis and min_ratio must be non-negative')
    gene, guide = data['gene'], data['guide']
    if not gene.obs_names.is_unique or not gene.obs_names.equals(guide.obs_names):
        raise ValueError('gene/guide cells must be unique and identically ordered')
    if not data.obs_names.equals(gene.obs_names):
        raise ValueError('top-level and RNA cell identities must agree')
    if 'guide_assignment' not in guide.layers:
        raise ValueError('guide_assignment layer is required; raw guide UMIs are not assignments')
    calls = sparse.csr_matrix(guide.layers['guide_assignment'], copy=True)
    if not np.isfinite(calls.data).all() or not np.isin(calls.data, [0, 1]).all():
        raise ValueError('guide_assignment must contain binary finite calls')
    calls.sum_duplicates()
    calls.eliminate_zeros()
    counts = np.asarray((calls > 0).sum(axis=1)).ravel()
    single = counts == 1

    # Raw guide UMIs of the assigned guide and of the strongest other guide, per cell.
    umis = sparse.csr_matrix(guide.X)
    if umis.shape != calls.shape:
        raise ValueError('guide counts and guide_assignment must have the same shape')
    if umis.nnz and (not np.isfinite(umis.data).all() or (umis.data < 0).any()):
        raise ValueError('guide counts must be finite and non-negative')
    assigned_umis = np.asarray(umis.multiply(calls).sum(axis=1)).ravel()
    others = umis - umis.multiply(calls)
    max_other_umis = others.max(axis=1).toarray().ravel() if others.nnz else np.zeros(umis.shape[0])

    enough_umis = assigned_umis >= min_umis
    dominant = assigned_umis >= min_ratio * max_other_umis
    keep = single & enough_umis & dominant
    if not keep.any():
        raise ValueError('No cell has exactly one assigned guide that passes the UMI floor and the ratio to the next guide')
    selected = data[gene.obs_names[keep], :].copy()
    # Nothing is refitted, collapsed or normalized. Every guide feature remains
    # available when counting assignments, before cis pair/gene subsetting.
    ratio_kept = assigned_umis[keep] / np.maximum(max_other_umis[keep], 1)
    report = {
        'experiment': 'high_moi_assignment_single_guide_inference',
        'input_cells': int(len(counts)),
        'selected_cells': int(keep.sum()),
        'removed_unassigned_cells': int((counts == 0).sum()),
        'removed_multi_guide_cells': int((counts > 1).sum()),
        'removed_low_assigned_umis': int((single & ~enough_umis).sum()),
        'removed_competing_guide': int((single & enough_umis & ~dominant).sum()),
        'min_assigned_umis': float(min_umis),
        'min_umi_ratio_to_next_guide': float(min_ratio),
        'selected_median_assigned_umis': float(np.median(assigned_umis[keep])),
        'selected_median_ratio_to_next_guide': float(np.median(ratio_kept)),
        'guide_features': int(guide.n_vars),
        'assignment_moi': str(guide.uns.get('moi', 'not recorded')),
        'inference_moi': 'low',
        'analysis_gene_rna_umis': int(selected['gene'].X.sum(dtype=np.int64)),
        'guide_counts_sha256': signature(selected['guide'].X),
        'guide_assignment_sha256': signature(selected['guide'].layers['guide_assignment']),
        'guide_counts_and_assignments_unchanged': True,
        'rna_sampling': False,
        'guide_reassignment': False,
        'source_qc': 'Inherited parent final eligibility; select rows by exactly one assigned guide that carries the UMI floor '
                     'and the UMI ratio to the next guide. RNA and guide counts/calls unchanged; no new RNA QC, doublet or HTO calling.',
    }
    return selected, report


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--output', default='single_guide_inference_input.h5mu', type=Path)
    p.add_argument('--report', default='inference_single_guide_filter.json', type=Path)
    p.add_argument('--cells', default='inference_single_guide_cells.tsv', type=Path)
    p.add_argument('--dataset', default='Inference cohort')
    p.add_argument('--min-umis', default=5, type=float, help='minimum raw UMIs of the assigned guide (default 5)')
    p.add_argument('--min-ratio', default=3.0, type=float,
                   help='minimum ratio of the assigned guide\'s UMIs to the largest other guide\'s UMIs (default 3)')
    args = p.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('Output must not overwrite the input')
    data = read_mudata_without_uns(args.input)
    # Avoid loading old multi-gigabyte inference tables on a resumed analysis,
    # but preserve cell annotations and all unrelated top-level metadata.
    with h5py.File(args.input, 'r') as f:
        if 'obs' in f:
            obs = read_elem(f['obs'])
            if not obs.index.equals(data.obs_names):
                raise ValueError('Stored top-level and modality cell identities differ')
            data.obs = obs
        for key, node in f.get('uns', {}).items():
            if key != 'pairs_to_test' and not key.endswith('_results'):
                data.uns[key] = read_elem(node)
    selected, report = select_single_guide_cells(data, min_umis=args.min_umis, min_ratio=args.min_ratio)
    report.update(dataset=args.dataset, source=str(args.input.resolve()))
    selected.write_h5mu(args.output)
    selected['gene'].obs_names.to_series(name='cell_id').to_csv(args.cells, sep='\t', index=False)
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
