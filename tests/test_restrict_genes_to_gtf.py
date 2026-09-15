"""Restricting the gene modality to the GTF's genes -- the TAP-seq panel filter.

Whole-transcriptome mapping of a targeted assay yields thousands of near-empty
off-panel genes; no expression threshold separates them from lowly expressed
panel genes (measured on the TAP-seq chr8 screen: 1% of cells drops 14 panel
genes and keeps 14 off-panel ones). The GTF is the panel, so the GTF is the
filter. These tests pin the id matching and the guard rails.
"""

import pathlib
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from create_mdata import restrict_genes_to_gtf  # noqa: E402


def _rna(var_names):
    n_cells, n_genes = 5, len(var_names)
    X = sparse.csr_matrix(np.arange(n_cells * n_genes, dtype=np.float32).reshape(n_cells, n_genes))
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_cells)]),
                       var=pd.DataFrame(index=var_names))
    adata.obs["percent_mito"] = np.linspace(0, 1, n_cells)
    return adata


def test_keeps_only_gtf_genes_and_all_cells():
    adata = _rna(["ENSG1", "ENSG2", "ENSG3", "ENSG4"])
    out = restrict_genes_to_gtf(adata, ["ENSG2", "ENSG4"])
    assert list(out.var_names) == ["ENSG2", "ENSG4"]
    assert out.n_obs == 5
    pd.testing.assert_series_equal(out.obs["percent_mito"], adata.obs["percent_mito"])
    # the kept columns are the original columns, not re-indexed copies
    np.testing.assert_array_equal(out.X.toarray(), adata[:, ["ENSG2", "ENSG4"]].X.toarray())


@pytest.mark.parametrize(
    "data_ids,gtf_ids",
    [
        (["ENSG1.5", "ENSG2.2", "ENSG3"], ["ENSG2", "ENSG3"]),      # versioned data, plain GTF
        (["ENSG1", "ENSG2", "ENSG3"], ["ENSG2.7", "ENSG3.1"]),      # plain data, versioned GTF
        (["ENSG1.1", "ENSG2.1", "ENSG3.1"], ["ENSG2.9", "ENSG3.3"]),  # both, different versions
    ],
)
def test_version_suffixes_are_ignored_on_both_sides(data_ids, gtf_ids):
    out = restrict_genes_to_gtf(_rna(data_ids), gtf_ids)
    assert [v.split(".")[0] for v in out.var_names] == ["ENSG2", "ENSG3"]


def test_gtf_genes_absent_from_the_data_are_tolerated(capsys):
    out = restrict_genes_to_gtf(_rna(["ENSG1", "ENSG2"]), ["ENSG2", "ENSG_NOT_MAPPED", "ENSG_ALSO_NOT"])
    assert list(out.var_names) == ["ENSG2"]
    assert "2 GTF genes absent from the data" in capsys.readouterr().out


def test_zero_overlap_is_an_error_not_an_empty_screen():
    with pytest.raises(ValueError, match="kept no genes"):
        restrict_genes_to_gtf(_rna(["ENSG1", "ENSG2"]), ["ENSGX", "ENSGY"])


def test_default_is_off_and_plumbed_through_both_call_sites():
    for name in ("nextflow.config", "nextflow_cc.config"):
        assert re.search(r"^\s*REFERENCE_restrict_genes_to_gtf\s*=\s*false", (REPO_ROOT / name).read_text(), re.M), name
    module = (REPO_ROOT / "modules/local/CreateMuData/main.nf").read_text()
    assert "val restrict_genes_to_gtf" in module and "--restrict-genes-to-gtf" in module
    workflow = (REPO_ROOT / "workflows/crispr_pipeline/main.nf").read_text()
    assert workflow.count("params.REFERENCE_restrict_genes_to_gtf,") == 2, "both CreateMuData call sites must pass it"
    # positional order in the module: capture_method, restrict flag, hashing file
    assert re.search(r"val capture_method\s+val restrict_genes_to_gtf\s+path adata_hashing", module)
