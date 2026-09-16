"""The pipeline's single cell filter, and the QC metric it must not destroy.

``QC_require_assigned_guide`` drops cells with no assigned guide in
``bin/mudata_concat.py``, so that both inference methods, every derived
covariate and both CRT pools see one cell population. The hard constraint is
the guide-assignment rate: recounted from the filtered object it reads 100% by
construction, so the counts as they stood before the filter are recorded in
``.uns`` and the QC path reports those. These tests pin the filter's boundary
(one guide is enough, none is not), the recorded counts, the off switch, the
JSON fields, and the fact that ``prepare_inference`` no longer narrows the cell
set a second time.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

anndata = pytest.importorskip("anndata")
mudata = pytest.importorskip("mudata")

import mudata_concat as mc  # noqa: E402
from qc_mudata_io import read_uns_scalars  # noqa: E402


# cell0: two guides, cell1: none, cell2: exactly one, cell3: none.
ASSIGNMENT = np.array(
    [
        [1, 1, 0],
        [0, 0, 0],
        [0, 0, 1],
        [0, 0, 0],
    ],
    dtype=np.uint16,
)
CELL_IDS = ["cell0", "cell1", "cell2", "cell3"]
GUIDE_IDS = ["g1", "g2", "nt1"]


def _mudata(assignment=ASSIGNMENT, with_layer=True, guide_umis=None):
    """A four-cell MuData whose guide calls are ASSIGNMENT."""
    guide_var = pd.DataFrame(
        {
            "guide_id": GUIDE_IDS,
            "intended_target_name": ["elem1", "elem1", "non-targeting"],
            "targeting": [True, True, False],
        },
        index=GUIDE_IDS,
    )
    guide_obs = pd.DataFrame(index=CELL_IDS)
    if guide_umis is not None:
        guide_obs["total_guide_umis"] = guide_umis
    guide = anndata.AnnData(
        # Raw guide UMIs are deliberately nonzero everywhere, so a test that
        # passes while reading guide.X instead of the assignment layer is not
        # passing by accident.
        X=sparse.csr_matrix(np.full_like(assignment, 3)),
        obs=guide_obs,
        var=guide_var,
    )
    if with_layer:
        guide.layers["guide_assignment"] = sparse.csr_matrix(assignment)

    gene = anndata.AnnData(
        X=sparse.csr_matrix(np.ones((len(CELL_IDS), 2), dtype=np.uint16)),
        obs=pd.DataFrame({"percent_mito": [1.0, 2.0, 3.0, 4.0]}, index=CELL_IDS),
        var=pd.DataFrame(index=["GENE1", "GENE2"]),
    )
    return mudata.MuData({"gene": gene, "guide": guide})


def test_keeps_a_cell_with_exactly_one_guide_and_drops_cells_with_none():
    filtered = mc.filter_cells_without_assigned_guide(_mudata())

    assert list(filtered.obs_names) == ["cell0", "cell2"]
    assert filtered["gene"].n_obs == 2
    assert filtered["guide"].n_obs == 2


def test_recorded_counts_match_the_pre_filter_reality():
    filtered = mc.filter_cells_without_assigned_guide(_mudata())

    assert filtered.uns["n_cells_before_assigned_guide_filter"] == 4
    assert filtered.uns["n_cells_after_assigned_guide_filter"] == 2
    assert filtered.uns["n_cells_with_assigned_guide"] == 2
    assert filtered.uns["n_cells_without_assigned_guide"] == 2
    assert filtered.uns["frac_cells_with_assigned_guide"] == 0.5
    assert filtered.uns["assigned_guide_filter_applied"] is True
    # The counts describe the population before the filter, not the object that
    # comes back out of it -- that is the whole point of recording them.
    assert filtered.n_obs == 2


def test_the_parameter_disables_the_filter_but_not_the_counts():
    kept = mc.filter_cells_without_assigned_guide(_mudata(), require_assigned_guide=False)

    assert list(kept.obs_names) == CELL_IDS
    assert kept.uns["assigned_guide_filter_applied"] is False
    assert kept.uns["n_cells_before_assigned_guide_filter"] == 4
    assert kept.uns["n_cells_after_assigned_guide_filter"] == 4
    assert kept.uns["n_cells_with_assigned_guide"] == 2
    assert kept.uns["n_cells_without_assigned_guide"] == 2


def test_counts_come_from_the_assignment_layer_not_raw_guide_umis():
    counts, source = mc.assigned_guides_per_cell(_mudata())

    assert list(counts) == [2, 0, 1, 0]
    assert source == "guide.layers['guide_assignment']"


def test_falls_back_to_guide_x_and_says_so(capsys):
    counts, source = mc.assigned_guides_per_cell(_mudata(with_layer=False))

    assert source == "guide.X"
    assert list(counts) == [3, 3, 3, 3]  # every raw UMI entry is nonzero
    assert "no 'guide_assignment' layer" in capsys.readouterr().out


def test_an_empty_result_is_refused_rather_than_written():
    empty = _mudata(assignment=np.zeros_like(ASSIGNMENT))
    with pytest.raises(ValueError, match="No cell carries an assigned guide"):
        mc.filter_cells_without_assigned_guide(empty)


def test_cells_are_filtered_before_genes(tmp_path):
    """The gene threshold is a fraction of the cells that are actually analysed."""
    mdata = _mudata()
    # GENE2 is detected only in the two cells that carry no guide.
    dense = np.array([[5, 0], [5, 7], [5, 0], [5, 7]], dtype=np.uint16)
    mdata.mod["gene"].X = sparse.csr_matrix(dense)

    filtered = mc.apply_cell_and_gene_filters(mdata, 0.5, require_assigned_guide=True)

    assert filtered.n_obs == 2
    # Over the retained cells GENE2 has zero support, so it goes; had the gene
    # filter run first it would have been detected in half the cells and stayed.
    assert list(filtered["gene"].var_names) == ["GENE1"]


def test_recorded_counts_survive_a_write_and_are_readable_without_uns(tmp_path):
    filtered = mc.filter_cells_without_assigned_guide(_mudata())
    path = tmp_path / "concat_mudata.h5mu"
    filtered.write(path)

    recorded = read_uns_scalars(path, mc.ASSIGNED_GUIDE_FILTER_KEYS)

    assert recorded["n_cells_before_assigned_guide_filter"] == 4
    assert recorded["n_cells_without_assigned_guide"] == 2
    assert recorded["frac_cells_with_assigned_guide"] == 0.5
    assert recorded["assigned_guide_filter_applied"] is True
    assert recorded["assigned_guide_filter_source"] == "guide.layers['guide_assignment']"


def test_guide_qc_reports_the_recorded_rate_not_a_recount():
    """Without this the assignment rate reads 100% on every filtered run."""
    from mapping_guide import compute_guide_metrics

    # The object the QC script sees: only the cells that survived the filter.
    survivors = mc.filter_cells_without_assigned_guide(_mudata())["guide"]
    survivors.obs["n_guides_per_cell"] = np.asarray(
        survivors.layers["guide_assignment"].sum(axis=1)
    ).ravel()

    recounted = compute_guide_metrics(survivors, include_per_guide_stats=False)
    assert recounted["frac_cells_with_guide"] == 1.0
    assert recounted["assigned_guide_counts_source"] == "final_mudata"

    reported = compute_guide_metrics(
        survivors,
        include_per_guide_stats=False,
        recorded_assignment_counts={
            "n_cells_before_assigned_guide_filter": 4,
            "n_cells_with_assigned_guide": 2,
            "n_cells_without_assigned_guide": 2,
            "frac_cells_with_assigned_guide": 0.5,
            "assigned_guide_filter_applied": True,
        },
    )
    assert reported["frac_cells_with_guide"] == 0.5
    assert reported["n_cells_with_guide"] == 2
    assert reported["n_cells_before_assigned_guide_filter"] == 4
    assert reported["n_cells_without_assigned_guide"] == 2
    assert reported["assigned_guide_counts_source"] == "mudata_concat"


def test_qc_metrics_json_carries_the_unassigned_count_and_fraction():
    from qc_metrics_json import _assigned_guide_filter_summary

    filtered = mc.filter_cells_without_assigned_guide(_mudata())
    summary = _assigned_guide_filter_summary(
        filtered, {"QC_require_assigned_guide": True}
    )

    assert summary["QC_require_assigned_guide"] is True
    assert summary["applied"] is True
    assert summary["cells_before"] == 4
    assert summary["cells_after"] == 2
    assert summary["cells_without_assigned_guide"] == 2
    assert summary["frac_cells_with_assigned_guide"] == 0.5


def test_qc_metrics_json_tolerates_a_mudata_written_before_the_filter():
    from qc_metrics_json import _assigned_guide_filter_summary

    summary = _assigned_guide_filter_summary(_mudata(), {})

    assert summary["cells_before"] is None
    assert summary["cells_without_assigned_guide"] is None
    assert "note" in summary


def test_the_derived_guide_umi_covariate_has_one_definition():
    """The adapter must consume the upstream column, not recompute its own."""
    from inference_covariates import (
        GUIDE_UMI_COVARIATE,
        derive_guide_umi_covariate,
        materialize_shared_covariates,
    )
    import perturbo_v2_pipeline_adapter as adapter

    prepared = mc.filter_cells_without_assigned_guide(
        _mudata(guide_umis=[10, 0, 40, 0])
    )
    materialize_shared_covariates(prepared)
    upstream = np.asarray(prepared["gene"].obs[GUIDE_UMI_COVARIATE], dtype=float)

    # Centred over the analysed cells only.
    expected = np.log1p([10.0, 40.0])
    np.testing.assert_allclose(upstream, expected - expected.mean())

    # The adapter leaves the upstream values alone ...
    adapter._ensure_covariates(prepared)
    np.testing.assert_allclose(
        np.asarray(prepared["gene"].obs[GUIDE_UMI_COVARIATE], dtype=float), upstream
    )

    # ... and derives exactly the same numbers when it has to compute them.
    fresh = mc.filter_cells_without_assigned_guide(_mudata(guide_umis=[10, 0, 40, 0]))
    adapter._ensure_covariates(fresh)
    np.testing.assert_allclose(
        np.asarray(fresh["gene"].obs[GUIDE_UMI_COVARIATE], dtype=float), upstream
    )
    assert derive_guide_umi_covariate(fresh) == GUIDE_UMI_COVARIATE


def test_the_derived_covariate_is_not_one_either_method_conditions_on():
    """It is written for provenance and the adapter's guard, not as a covariate."""
    from inference_covariates import (
        CANONICAL_COVARIATES,
        GUIDE_UMI_COVARIATE,
        materialize_shared_covariates,
    )

    assert GUIDE_UMI_COVARIATE not in {name for name, _kind, _aliases in CANONICAL_COVARIATES}

    prepared = mc.filter_cells_without_assigned_guide(_mudata(guide_umis=[10, 0, 40, 0]))
    written = materialize_shared_covariates(prepared)

    # SCEPTRE reads the top-level obs as colData; the guide-UMI term is not there.
    assert GUIDE_UMI_COVARIATE not in written
    assert GUIDE_UMI_COVARIATE not in prepared.obs.columns


def _cis_subset_mudata():
    """A four-cell screen where one cell carries only an untested, non-control guide.

    That cell is the regression: ``prepare_inference`` used to drop it, because
    it fell out of the cis guide subset, which is how SCEPTRE ended up analysing
    a different cell set from PerTurbo.
    """
    guide_ids = ["g1", "nt1", "g_untested"]
    guide_var = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "intended_target_name": ["ELEM1", "non-targeting", "ELEM9"],
            "intended_target_chr": ["chr1", "", "chr9"],
            "intended_target_start": [100.0, np.nan, 900.0],
            "intended_target_end": [150.0, np.nan, 950.0],
            "targeting": [True, False, True],
            "type": ["targeting", "non-targeting", "targeting"],
        },
        index=guide_ids,
    )
    assignment = np.array(
        [
            [1, 0, 0],  # cell0: the tested guide
            [0, 1, 0],  # cell1: a control guide
            [0, 0, 1],  # cell2: only an untested, non-control guide
            [1, 0, 0],  # cell3: the tested guide
        ],
        dtype=np.uint16,
    )
    guide = anndata.AnnData(
        X=sparse.csr_matrix(assignment),
        obs=pd.DataFrame(index=CELL_IDS),
        var=guide_var,
    )
    guide.layers["guide_assignment"] = sparse.csr_matrix(assignment)

    gene = anndata.AnnData(
        X=sparse.csr_matrix(np.ones((len(CELL_IDS), 2), dtype=np.uint16)),
        obs=pd.DataFrame({"percent_mito": [1.0, 2.0, 3.0, 4.0]}, index=CELL_IDS),
        var=pd.DataFrame(index=["GENE1", "GENE2"]),
    )
    return mudata.MuData({"gene": gene, "guide": guide})


def test_prepare_inference_subsets_genes_and_guides_but_not_cells(tmp_path, monkeypatch):
    import prepare_inference

    mudata_path = tmp_path / "concat_mudata.h5mu"
    _cis_subset_mudata().write(mudata_path)

    pairs_path = tmp_path / "pairs_to_test.csv"
    pd.DataFrame({"guide_id": ["g1"], "gene_name": ["GENE1"]}).to_csv(
        pairs_path, index=False
    )

    monkeypatch.chdir(tmp_path)
    prepare_inference.main(str(pairs_path), str(mudata_path), subset_for_cis=True)

    prepared = mudata.read_h5mu(tmp_path / "mudata_inference_input.h5mu")

    # Genes and guides are still narrowed to the cis set plus controls ...
    assert list(prepared["gene"].var_names) == ["GENE1"]
    assert sorted(prepared["guide"].var["guide_id"].astype(str)) == ["g1", "nt1"]

    # ... and the cell set is untouched, including the cell whose only guide was
    # dropped by the guide subset.
    assert list(prepared.obs_names) == CELL_IDS
    assert prepared["gene"].n_obs == 4
    assert prepared["guide"].n_obs == 4

    # The shared covariates are written here, over those same cells, and the
    # derived guide-UMI column survives the write for the adapter to consume.
    assert "percent_mito" in prepared.obs.columns
    assert "log1p_total_guide_umis_centered" in prepared["gene"].obs.columns


def test_prepare_inference_no_longer_mentions_a_cell_filter():
    source = (BIN_DIR / "prepare_inference.py").read_text()
    assert "targeted_cells" not in source
