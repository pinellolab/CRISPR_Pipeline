import importlib.util
from pathlib import Path
import sys

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
spec = importlib.util.spec_from_file_location('postconcat_embedding_qc', ROOT / 'bin/postconcat_embedding_qc.py')
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


def fixture_mudata():
    cells = [f'cell_{i}' for i in range(1000)]
    genes = ['ENSG_GLOBAL.1', 'ENSG_RESCUE.2', 'ENSG_TOO_RARE.1', 'ENSG_CONTROL.1', 'ENSG_ZERO.1']
    matrix = sp.lil_matrix((1000, len(genes)), dtype=np.uint16)
    matrix[:50, 0] = 1
    matrix[:10, 1] = 1
    matrix[0, 2] = 1
    matrix[:10, 3] = 1
    gene = ad.AnnData(
        matrix.tocsr(),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(
            {'symbol': ['GLOBAL', 'RESCUE', 'RARE', 'CONTROL', 'ZERO']},
            index=genes,
        ),
    )

    guide_var = pd.DataFrame({
        'guide_id': ['g_rescue_a', 'g_rescue_b', 'g_rare', 'g_control', 'g_zero'],
        'targeting': [True, True, True, True, True],
        'type': ['targeting', 'targeting', 'targeting', 'negative control', 'targeting'],
        'intended_target_name': [
            'ENSG_RESCUE', 'ENSG_RESCUE', 'ENSG_TOO_RARE', 'ENSG_CONTROL', 'ENSG_ZERO'
        ],
    }, index=['g_rescue_a', 'g_rescue_b', 'g_rare', 'g_control', 'g_zero'])
    assignment = sp.lil_matrix((1000, 5), dtype=np.uint8)
    assignment[:12, 0] = 1
    assignment[8:20, 1] = 1
    assignment[:20, 2] = 1
    assignment[:100, 3] = 1
    assignment[:100, 4] = 1
    guide = ad.AnnData(
        sp.csr_matrix((1000, 5)),
        obs=pd.DataFrame(index=cells),
        var=guide_var,
    )
    guide.layers['guide_assignment'] = assignment.tocsr()
    return md.MuData({'gene': gene, 'guide': guide})


def test_target_aware_rule_is_inclusive_and_excludes_controls_and_zero_counts():
    mdata = fixture_mudata()
    keep, required, table, details = qc.gene_filter_decisions(
        mdata, mdata.mod['gene'], 0.05, True, 20, 0.001,
    )
    assert required == 50
    assert keep.tolist() == [True, True, True, False, False]
    assert table['assigned_target_positive_cells'].tolist() == [0, 20, 20, 0, 100]
    assert table['keep_reason'].tolist() == [
        'global_prevalence', 'target_aware_rescue', 'target_aware_rescue', 'removed', 'removed'
    ]
    assert details['canonical_intended_target_genes'] == 3


def test_target_rescue_is_opt_in():
    mdata = fixture_mudata()
    keep, required, table, details = qc.gene_filter_decisions(
        mdata, mdata.mod['gene'], 0.05, False, 20, 0.001,
    )
    assert required == 50
    assert keep.tolist() == [True, False, False, False, False]
    assert not details
