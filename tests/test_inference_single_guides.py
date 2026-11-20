"""Inference row selection must use assigned guides, apply the UMI floor and the ratio to the next guide, and preserve raw evidence."""
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from filter_inference_single_guides import select_single_guide_cells


def fixture():
    # zero: no call; target: one call, 10 UMIs, next guide 3 (ratio 3.3); multi: two calls; ntc: one call, 7 UMIs, no other guide;
    # faint: one call at 4 UMIs (under the floor); crowded: one call at 12 UMIs with a 6-UMI competitor (ratio 2);
    # background: one call at 20 UMIs over 1-UMI background guides (ratio 20).
    obs = pd.DataFrame(index=['zero', 'target', 'multi', 'ntc', 'faint', 'crowded', 'background'])
    gene = ad.AnnData(sparse.csr_matrix([[8, 2], [5, 9], [7, 4], [9, 2], [6, 6], [3, 8], [4, 4]]), obs=obs.copy())
    guide = ad.AnnData(sparse.csr_matrix([[20, 8, 0], [10, 3, 0], [5, 9, 0], [0, 0, 7], [4, 0, 1], [12, 6, 0], [20, 1, 1]]), obs=obs.copy())
    guide.layers['guide_assignment'] = sparse.csr_matrix([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0], [1, 0, 0]])
    guide.uns['moi'] = 'high'
    return md.MuData({'gene': gene, 'guide': guide})


def test_assigned_not_raw_counts_and_no_reassignment():
    data = fixture()
    got, report = select_single_guide_cells(data)
    assert list(got.obs_names) == ['target', 'ntc', 'background']
    assert data.n_obs == 7 and got['guide'].n_vars == 3
    assert got['guide'].uns['moi'] == 'high' and report['inference_moi'] == 'low'
    for mod in ['gene', 'guide']:
        assert (got[mod].X != data[mod][['target', 'ntc', 'background']].X).nnz == 0
    assert got['guide'].layers['guide_assignment'].toarray().tolist() == [[1, 0, 0], [0, 0, 1], [1, 0, 0]]
    assert report['removed_unassigned_cells'] == report['removed_multi_guide_cells'] == 1
    assert report['removed_low_assigned_umis'] == 1 and report['removed_competing_guide'] == 1
    assert report['min_assigned_umis'] == 5 and report['min_umi_ratio_to_next_guide'] == 3
    assert report['selected_median_assigned_umis'] == 10
    assert report['selected_cells'] == 3 and report['input_cells'] == 7


def test_thresholds_are_the_only_difference_from_exactly_one_call():
    got, report = select_single_guide_cells(fixture(), min_umis=0, min_ratio=0)
    assert list(got.obs_names) == ['target', 'ntc', 'faint', 'crowded', 'background']
    assert report['removed_low_assigned_umis'] == report['removed_competing_guide'] == 0
    got, _ = select_single_guide_cells(fixture(), min_umis=0, min_ratio=3)
    assert list(got.obs_names) == ['target', 'ntc', 'faint', 'background']
    got, _ = select_single_guide_cells(fixture(), min_umis=5, min_ratio=0)
    assert list(got.obs_names) == ['target', 'ntc', 'crowded', 'background']


def test_ratio_requires_the_assigned_guide_to_lead():
    data = fixture()
    # the call sits on a 3-UMI guide while a 10-UMI guide is uncalled: a mis-called cell, rejected by the ratio
    data['guide'].layers['guide_assignment'] = sparse.csr_matrix([[0, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0], [1, 0, 0]])
    got, report = select_single_guide_cells(data, min_umis=0, min_ratio=1)
    assert 'target' not in got.obs_names and report['removed_competing_guide'] == 1


def test_fail_closed():
    data = fixture(); del data['guide'].layers['guide_assignment']
    with pytest.raises(ValueError, match='raw guide UMIs'):
        select_single_guide_cells(data)
    data = fixture(); data['guide'].layers['guide_assignment'] = sparse.csr_matrix(np.zeros((7, 3)))
    with pytest.raises(ValueError, match='No cell'):
        select_single_guide_cells(data)
    data = fixture(); data['guide'].layers['guide_assignment'] = data['guide'].layers['guide_assignment'].astype(float)
    data['guide'].layers['guide_assignment'][1, 0] = 0.8
    with pytest.raises(ValueError, match='binary'):
        select_single_guide_cells(data)
    with pytest.raises(ValueError, match='non-negative'):
        select_single_guide_cells(fixture(), min_ratio=-1)


def test_round_trip_preserves_counts_and_unrelated_metadata(tmp_path):
    data = fixture(); data.uns['unrelated_metadata'] = 'keep me'
    data.uns['global_analysis_per_guide_results'] = pd.DataFrame({'old': [1]})
    source = tmp_path / 'parent.h5mu'; data.write_h5mu(source)
    subprocess.run([sys.executable, str(ROOT / 'bin/filter_inference_single_guides.py'), '--input', str(source), '--min-umis', '5', '--min-ratio', '3'],
                   cwd=tmp_path, check=True)
    got = md.read_h5mu(tmp_path / 'single_guide_inference_input.h5mu')
    assert got.uns['unrelated_metadata'] == 'keep me'
    assert 'global_analysis_per_guide_results' not in got.uns
    assert list(got.obs_names) == ['target', 'ntc', 'background']
    assert (got['gene'].X != data['gene'][['target', 'ntc', 'background']].X).nnz == 0
    report = json.loads((tmp_path / 'inference_single_guide_filter.json').read_text())
    assert report['selected_cells'] == 3 and report['guide_counts_and_assignments_unchanged'] is True
    cells = pd.read_csv(tmp_path / 'inference_single_guide_cells.tsv', sep='\t')
    assert cells['cell_id'].tolist() == ['target', 'ntc', 'background']
