import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.mark.parametrize('shared,extra,legacy,expected', [
    (30, 300, 0, 'PASS'), (10, 300, 0, 'FAIL'),
    (0, 300, 0, 'FAIL'), (30, 300, 0.5, 'FAIL'),
])
def test_filtered_rna_overlap(tmp_path, monkeypatch, shared, extra, legacy, expected):
    monkeypatch.setenv('MPLCONFIGDIR', str(tmp_path / 'matplotlib'))
    script = Path(__file__).resolve().parents[1] / 'bin/guide_mapping_qc.py'
    spec = importlib.util.spec_from_file_location('guide_mapping_qc', script)
    qc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(qc)
    def matrix(names):
        return ad.AnnData(sparse.csr_matrix(np.ones((len(names), 2))),
                          obs=pd.DataFrame({'batch': ['set1'] * len(names)}, index=names),
                          var=pd.DataFrame(index=['g1', 'g2']))
    rna = matrix([f'cell{i}' for i in range(30)])
    guide = matrix([f'cell{i}' for i in range(shared)] + [f'raw{i}' for i in range(extra)])
    rna.write_h5ad(tmp_path / 'rna.h5ad')
    guide.write_h5ad(tmp_path / 'guide.h5ad')
    md.MuData({'gene': rna}).write_h5mu(tmp_path / 'data.h5mu')
    (tmp_path / 'guides.tsv').write_text('guide_id\ng1\ng2\n')
    monkeypatch.setattr(sys, 'argv', [str(script), '--rna', str(tmp_path / 'rna.h5ad'),
        '--guide', str(tmp_path / 'guide.h5ad'), '--mudata', str(tmp_path / 'data.h5mu'),
        '--guide-metadata', str(tmp_path / 'guides.tsv'), '--reverse-complement-guides', 'true',
        '--min-overlap-to-guide-fraction', str(legacy), '--outdir', str(tmp_path / 'qc')])
    qc.main()
    report = json.loads((tmp_path / 'qc/guide_mapping_qc.json').read_text())
    assert report['status'] == expected
    assert report['measurement_sets'][0]['overlap_to_rna_fraction'] == shared / 30
    assert report['overall']['recovered_guide_fraction'] == 1
