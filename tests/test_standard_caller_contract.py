"""Keep native/ambient caller experiments outside the standard run contract."""
import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from inference_covariates import validate_positive


def test_preprocessing_exposes_only_standard_barcode_selectors():
    tree = ast.parse((ROOT / 'bin/preprocess_adata.py').read_text())
    choices = [
        ast.literal_eval(keyword.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and node.args
        and isinstance(node.args[0], ast.Constant)
        and node.args[0].value == '--barcode-filter'
        for keyword in node.keywords if keyword.arg == 'choices'
    ]
    assert choices == [['none', 'knee', 'knee2']]
    for name in ['nextflow.config', 'nextflow_schema.json',
                 'modules/local/PreprocessAnnData/main.nf',
                 'subworkflows/local/preprocessing_pipeline/main.nf',
                 'bin/preprocess_adata.py']:
        source = (ROOT / name).read_text().lower()
        for experimental in ['cell_ranger_empty_drops', 'cell_calling_chemistry',
                             'native_cell_calls', 'rna_emptydrops', 'ordmag']:
            assert experimental not in source, (name, experimental)


def test_zero_guide_covariate_is_finite_without_changing_input_counts():
    raw = np.array([0., 1., 7.])
    transformed = validate_positive(raw, 'total_guide_umis')
    np.testing.assert_array_equal(raw, [0., 1., 7.])
    np.testing.assert_array_equal(transformed, [1., 1., 7.])
    assert np.isfinite(np.log(transformed)).all()


@pytest.mark.parametrize('bad', [-1., np.nan, np.inf])
def test_invalid_guide_counts_are_rejected(bad):
    with pytest.raises(ValueError):
        validate_positive([bad], 'total_guide_umis')


@pytest.mark.parametrize('column', ['total_gene_umis', 'num_expressed_genes'])
def test_rna_zero_depth_stays_an_error(column):
    with pytest.raises(ValueError):
        validate_positive([0.], column)
