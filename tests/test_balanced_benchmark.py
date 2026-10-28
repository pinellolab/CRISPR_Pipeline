"""Scientific invariants of matching raw RNA while freezing guide evidence."""
import importlib.util
import json
import sys
import subprocess
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
spec = importlib.util.spec_from_file_location('balance_mudata', ROOT / 'bin/balance_mudata.py')
balance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(balance)


def test_exact_joint_hypergeometric():
    raw = sparse.csr_matrix([[0, 12, 8], [9, 0, 3], [1, 7, 2]])
    a, total = balance.exact_thin(raw, 24, 42)
    b, _ = balance.exact_thin(raw, 24, 42)
    assert total == 42
    assert a.sum() == 24
    assert (a != b).nnz == 0
    assert np.all(a.toarray() <= raw.toarray())
    assert (raw.toarray() == [[0, 12, 8], [9, 0, 3], [1, 7, 2]]).all()
    for invalid in (sparse.csr_matrix([[-1, 3]]), sparse.csr_matrix([[0.5, 3]])):
        with pytest.raises(ValueError): balance.exact_thin(invalid, 1, 42)
    with pytest.raises(ValueError): balance.exact_thin(raw, 43, 42)
    selected = balance.select_cells(['c', 'b', 'a'], 2, 42)
    assert selected == balance.select_cells(['a', 'c', 'b'], 2, 42)


def test_real_mudata_round_trip(tmp_path, monkeypatch):
    cells = pd.Index([f'cell{i}' for i in range(12)])
    rna = ad.AnnData(sparse.csr_matrix(np.arange(48).reshape(12, 4) % 8 + 1),
                     obs=pd.DataFrame(index=cells),
                     var=pd.DataFrame({'gene_name': ['MT-A', 'RPL1', 'ABC', 'XYZ']}, index=['g0','g1','g2','g3']))
    rna.layers['normalized'] = rna.X.copy()
    rna.obsm['X_pca'] = np.zeros((12, 2))
    guide = ad.AnnData(sparse.csr_matrix(np.arange(24).reshape(12, 2) % 3),
                       obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=['q0', 'q1']))
    guide.layers['guide_assignment'] = sparse.csr_matrix(np.tile([[1, 0]], (12, 1)))
    data = md.MuData({'gene': rna, 'guide': guide})
    source = tmp_path/'source.h5mu'; parent = tmp_path/'parent.h5mu'
    data.write_h5mu(source)
    baseline = data[:10].copy()
    baseline.mod['gene'] = baseline['gene'][:, ['g0','g2','g3']].copy()
    baseline.update()
    baseline.uns['old_inference_results'] = {'pvalue': 0.01}
    baseline.write_h5mu(parent)
    out = tmp_path/'matched'
    monkeypatch.setattr(sys, 'argv', ['balance', '--source',str(source), '--parent-inference',str(parent),
                                     '--dataset','toy', '--cells','8', '--rna-umis','60','--seed','42','--outdir',str(out)])
    balance.main()
    full = md.read_h5mu(out/'sampled_full_mudata.h5mu')
    analysis = md.read_h5mu(out/'balanced_inference_input.h5mu')
    selected = balance.select_cells(cells[:10], 8, 42)
    assert full.n_obs == analysis.n_obs == 8
    assert list(full.obs_names) == selected
    assert full['gene'].X.sum() == 60
    assert list(analysis['gene'].var_names) == ['g0','g2','g3']
    assert not [key for key in full['gene'].layers if key is not None] and not full['gene'].obsm
    assert 'old_inference_results' not in analysis.uns
    expected = data['guide'][selected]
    assert (full['guide'].X != expected.X).nnz == 0
    assert (full['guide'].layers['guide_assignment'] != expected.layers['guide_assignment']).nnz == 0
    assert np.array_equal(full['gene'].obs.total_gene_umis, np.asarray(full['gene'].X.sum(1)).ravel())
    report = json.loads((out/'downsampling_report.json').read_text())
    assert report['observed_total_rna_umis'] == 60
    assert report['analysis_gene_rna_umis'] == analysis['gene'].X.sum()
    with pytest.raises(ValueError, match='overwrite'): balance.main()
    # Exercise the final integrity gate and scientific report without loading
    # or copying parent result tables into this synthetic fixture.
    dashboard = tmp_path/'dashboard'; dashboard.mkdir()
    for folder in ('additional_qc','evaluation_output','control_output','benchmark_output'):
        (dashboard/folder).mkdir()
    (dashboard/'params.json').write_text('{}')
    subprocess.run([sys.executable,str(ROOT/'bin/balanced_dashboard.py'),
                    '--input',str(out/'balanced_inference_input.h5mu'),
                    '--report',str(out/'downsampling_report.json'),
                    '--params',str(dashboard/'params.json')],cwd=dashboard,check=True)
    metrics=json.loads((dashboard/'pipeline_qc_metrics.json').read_text())
    assert metrics['balanced_matching']['observed_total_rna_umis'] == 60
    assert (dashboard/'dashboard.html').is_file()


def test_full_workflow_contract():
    flow = (ROOT/'workflows/balanced_benchmark/main.nf').read_text()
    assert 'inference_pipeline(Balanced.inference_input' in flow
    assert 'evaluation_pipeline' in flow and 'balanced_dashboard(' in flow
    assert 'Balanced.full_mudata' in flow and 'Inference.calibration_output' in flow
    assert 'INFERENCE_input_mudata =' not in flow


def test_normalization_preserves_raw_input(tmp_path):
    rng = np.random.default_rng(71)
    counts = sparse.csr_matrix(rng.poisson(1.5, (80, 100)))
    obs = pd.DataFrame({'batch': ['a']*40 + ['b']*40}, index=[f'c{i}' for i in range(80)])
    rna = ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=[f'g{i}' for i in range(100)]))
    balance.update_rna(rna)
    guide = ad.AnnData(sparse.csr_matrix(np.ones((80, 2))), obs=obs.copy())
    guide.layers['guide_assignment'] = guide.X.copy()
    source = tmp_path/'normalization_input.h5mu'
    md.MuData({'gene':rna, 'guide':guide}).write_h5mu(source)
    before = balance.sha256(source)
    subprocess.run([sys.executable, str(ROOT/'bin/balanced_normalization.py'), '--input',str(source),
                    '--outdir',str(tmp_path/'plots'), '--cell-cycle-genes',
                    str(ROOT/'assets/cell_cycle/regev_lab_cell_cycle_genes.txt')], check=True)
    report = json.loads((tmp_path/'plots/embedding_qc_metrics.json').read_text())
    assert balance.sha256(source) == before
    assert report['retained_cells'] == 80 and report['normalized_matrix_saved'] is False
    assert (tmp_path/'plots/normalization_check.png').is_file()
    assert (tmp_path/'plots/umap_qc_panel.png').is_file()
