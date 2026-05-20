import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / 'bin'
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import build_catalog_per_element_output as catalog_builder


def _make_test_mudata():
    obs = pd.DataFrame(index=['cell1', 'cell2', 'cell3'])

    gene_var = pd.DataFrame({'symbol': ['SYM1', 'SYM2']}, index=['GENE1', 'GENE2'])
    gene = ad.AnnData(X=np.array([[1, 0], [0, 1], [1, 1]], dtype=float), obs=obs.copy(), var=gene_var)

    guide_var = pd.DataFrame(
        {
            'guide_id': ['gA2', 'gA1', 'gB1'],
            'intended_target_name': ['elemA', 'elemA', 'elemB'],
            'intended_target_chr': ['chr1', 'chr1', 'chr2'],
            'intended_target_start': [100, 100, 300],
            'intended_target_end': [200, 200, 400],
            'type': ['targeting', 'targeting', 'targeting'],
        },
        index=['gA2', 'gA1', 'gB1'],
    )
    guide = ad.AnnData(
        X=np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ],
            dtype=float,
        ),
        obs=obs.copy(),
        var=guide_var,
    )
    guide.layers['guide_assignment'] = sparse.csr_matrix(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ],
            dtype=float,
        )
    )

    return mu.MuData({'gene': gene, 'guide': guide})


def _make_input_tables():
    cis = pd.DataFrame(
        {
            'gene_id': ['GENE1'],
            'intended_target_name': ['elemA'],
            'intended_target_chr': ['chr1'],
            'intended_target_start': [100],
            'intended_target_end': [200],
            'sceptre_log2_fc': [1.2],
            'sceptre_p_value': [0.01],
        }
    )

    trans = pd.DataFrame(
        {
            'gene_id': ['GENE1', 'GENE2'],
            'intended_target_name': ['elemA', 'elemB'],
            'intended_target_chr': ['chr1', 'chr2'],
            'intended_target_start': [100, 300],
            'intended_target_end': [200, 400],
            'log2_fc': [0.5, -1.0],
            'p_value': [0.2, 0.0],
        }
    )
    return cis, trans


def test_create_catalog_per_element_merge_and_metrics():
    mdata = _make_test_mudata()
    cis, trans = _make_input_tables()

    catalog = catalog_builder.create_catalog_per_element(cis, trans, mdata)

    assert list(catalog.columns) == catalog_builder.OUTPUT_COLUMNS
    assert len(catalog) == 2

    row_a = catalog[(catalog['element_name'] == 'elemA') & (catalog['gene_id'] == 'GENE1')].iloc[0]
    assert row_a['sceptre_log2_fc'] == 1.2
    assert np.isclose(row_a['sceptre_log10_p_value'], 2.0)
    assert row_a['perturbo_log2_fc'] == 0.5
    assert np.isclose(row_a['perturbo_log10_p_value'], -np.log10(0.2))
    assert row_a['guide_ids'] == 'gA1;gA2'
    assert row_a['gene_name'] == 'SYM1'
    assert row_a['nPerturbedCells'] == 2

    row_b = catalog[(catalog['element_name'] == 'elemB') & (catalog['gene_id'] == 'GENE2')].iloc[0]
    assert pd.isna(row_b['sceptre_log2_fc'])
    assert pd.isna(row_b['sceptre_log10_p_value'])
    assert row_b['perturbo_log2_fc'] == -1.0
    assert row_b['perturbo_log10_p_value'] == 300.0
    assert row_b['gene_name'] == 'SYM2'
    assert row_b['nPerturbedCells'] == 1


def test_build_catalog_writes_gzip_tsv(tmp_path):
    mdata = _make_test_mudata()
    cis, trans = _make_input_tables()

    cis_path = tmp_path / 'cis_per_element_output.tsv.gz'
    trans_path = tmp_path / 'trans_per_element_output.tsv.gz'
    mudata_path = tmp_path / 'inference_mudata.h5mu'
    out_path = tmp_path / 'catalog_per_element_output.tsv.gz'

    cis.to_csv(cis_path, sep='\t', index=False, compression='gzip')
    trans.to_csv(trans_path, sep='\t', index=False, compression='gzip')
    mdata.write(mudata_path)

    catalog_builder.build_catalog_per_element_output(
        cis_per_element_path=str(cis_path),
        trans_per_element_path=str(trans_path),
        mudata_path=str(mudata_path),
        output_path=str(out_path),
    )

    assert out_path.exists()
    observed = pd.read_csv(out_path, sep='\t')
    assert list(observed.columns) == catalog_builder.OUTPUT_COLUMNS
    assert not observed.empty
