import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import perturbo_v2_pipeline_adapter as adapter


def _make_mudata():
    obs = pd.DataFrame(index=["cell1", "cell2", "cell3"])
    gene = ad.AnnData(
        X=np.array([[10, 1], [2, 8], [4, 4]], dtype=np.float32),
        obs=obs.copy(),
        var=pd.DataFrame({"symbol": ["SYM1", "SYM2"]}, index=["GENE1", "GENE2"]),
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": ["gA", "gB", "nt1"],
            "targeting": [True, True, False],
            "type": ["targeting", "targeting", "non-targeting"],
            "intended_target_name": ["elemA", "elemB", "non-targeting"],
            "intended_target_chr": ["chr1", "chr2", ""],
            "intended_target_start": [100.0, 300.0, np.nan],
            "intended_target_end": [200.0, 400.0, np.nan],
        },
        index=["gA", "gB", "nt1"],
    )
    guide = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
                dtype=np.float32,
            )
        ),
        obs=obs.copy(),
        var=guide_var,
    )
    guide.layers["guide_assignment"] = guide.X.copy()
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["pairs_to_test"] = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2", "GENE1"],
            "guide_id": ["gA", "gB", "nt1"],
        }
    )
    return mdata


def test_prepare_mudata_adds_v2_metadata_and_control_guide_names(tmp_path):
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _make_mudata().write(input_path)

    guide_name_map = adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)

    prepared = mu.read_h5mu(prepared_path)
    guide = prepared["guide"]
    gene = prepared["gene"]

    assert adapter.ELEMENT_MAP_KEY in guide.varm
    assert adapter.ELEMENT_NAMES_KEY in guide.uns
    assert adapter.GUIDE_MAP_KEY in guide.varm
    assert adapter.GUIDE_NAMES_KEY in guide.uns
    assert adapter.ELEMENT_PAIR_MASK_KEY in gene.varm
    assert adapter.GUIDE_PAIR_MASK_KEY in gene.varm
    assert "total_gene_umis" in gene.obs
    assert "log1p_total_guide_umis_centered" in gene.obs
    assert any(str(name).startswith("non-targeting|") for name in guide.uns[adapter.GUIDE_NAMES_KEY])
    assert guide_name_map["non-targeting|nt1"] == "nt1"
    element_mask = sparse.csr_matrix(gene.varm[adapter.ELEMENT_PAIR_MASK_KEY])
    guide_mask = sparse.csr_matrix(gene.varm[adapter.GUIDE_PAIR_MASK_KEY])
    assert element_mask.shape[0] == gene.n_vars
    assert guide_mask.shape[0] == gene.n_vars
    assert element_mask.nnz >= 3
    assert guide_mask.nnz >= 3


def test_convert_element_effects_maps_metadata_and_filters_pairs(tmp_path):
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _make_mudata().write(input_path)
    adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)
    prepared = mu.read_h5mu(prepared_path)
    elem_a = prepared["guide"].var.loc["gA", "intended_target_key"]
    elem_b = prepared["guide"].var.loc["gB", "intended_target_key"]

    effects = pd.DataFrame(
        {
            "element": [elem_a, elem_a, elem_b],
            "gene": ["GENE1", "GENE2", "GENE2"],
            "posterior_mean": [np.log(2), np.log(4), -np.log(2)],
            "posterior_scale": [np.log(2) / 10, np.log(2) / 5, np.log(2) / 2],
            "posterior_prob": [0.3, 0.4, 0.5],
            "empirical_p_value": [0.01, np.nan, 0.2],
        }
    )

    observed = adapter.convert_element_effects(
        effects,
        prepared_path,
        test_all_pairs=False,
    )

    assert list(observed["gene_id"]) == ["GENE1", "GENE2"]
    assert list(observed["intended_target_name"]) == ["elemA", "elemB"]
    assert np.isclose(observed.loc[0, "log2_fc"], 1.0)
    assert np.isclose(observed.loc[0, "perturbo_fc_se"], 0.1)
    assert np.isclose(observed.loc[0, "p_value"], 0.01)
    assert "perturbo_q_value" in observed.columns


def test_convert_guide_effects_restores_control_guide_ids_and_filters_pairs(tmp_path):
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _make_mudata().write(input_path)
    guide_name_map = adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)

    effects = pd.DataFrame(
        {
            "element": ["gA", "gA", "non-targeting|nt1"],
            "gene": ["GENE1", "GENE2", "GENE1"],
            "posterior_mean": [np.log(2), np.log(3), 0.0],
            "posterior_scale": [np.log(2) / 10, np.log(2) / 9, np.log(2)],
            "posterior_prob": [0.05, 0.07, 0.8],
            "empirical_p_value": [np.nan, np.nan, 0.9],
        }
    )

    observed = adapter.convert_guide_effects(
        effects,
        guide_name_map,
        prepared_path,
        test_all_pairs=False,
    )

    assert list(observed["guide_id"]) == ["gA", "nt1"]
    assert list(observed["gene_id"]) == ["GENE1", "GENE1"]
    assert np.isclose(observed.loc[0, "log2_fc"], 1.0)
    assert np.isclose(observed.loc[0, "p_value"], 0.05)
    assert np.isclose(observed.loc[1, "p_value"], 0.9)


@pytest.mark.parametrize(
    ("test_all_pairs", "analysis_prefix"),
    [(False, "local_analysis"), (True, "global_analysis")],
)
def test_run_pipeline_adapter_patches_uns_without_rewriting_x(
    tmp_path, monkeypatch, test_all_pairs, analysis_prefix
):
    pytest.importorskip("pyarrow")
    input_path = tmp_path / "input.h5mu"
    _make_mudata().write(input_path)

    import h5py

    with h5py.File(input_path, "r") as f:
        x_before = f["mod/gene/X"][:].copy()

    # Element-level "element" values are intended_target_key groupings, not
    # raw guide ids -- derive them the same way convert_element_effects does.
    prepared_path = tmp_path / "prepared_for_test.h5mu"
    adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)
    prepared = mu.read_h5mu(prepared_path)
    elem_a = prepared["guide"].var.loc["gA", "intended_target_key"]
    elem_b = prepared["guide"].var.loc["gB", "intended_target_key"]

    element_effects = pd.DataFrame(
        {
            "element": [elem_a, elem_b],
            "gene": ["GENE1", "GENE2"],
            "posterior_mean": [np.log(2), np.log(3)],
            "posterior_scale": [np.log(2) / 10, np.log(2) / 9],
            "posterior_prob": [0.05, 0.07],
        }
    )
    guide_effects = pd.DataFrame(
        {
            "element": ["gA", "gB"],
            "gene": ["GENE1", "GENE2"],
            "posterior_mean": [np.log(2), np.log(3)],
            "posterior_scale": [np.log(2) / 10, np.log(2) / 9],
            "posterior_prob": [0.05, 0.07],
        }
    )

    def fake_run_perturbo(input_path, out_dir, *, map_key, names_key, args, **kwargs):
        out_dir.mkdir(parents=True, exist_ok=True)
        effects = element_effects if map_key == adapter.ELEMENT_MAP_KEY else guide_effects
        effects.to_parquet(out_dir / "element_effects.parquet", index=False)
        return object()

    monkeypatch.setattr(adapter, "_run_perturbo", fake_run_perturbo)
    monkeypatch.setattr(adapter, "_wait_for_fits", lambda processes: None)

    output_mudata = tmp_path / "output.h5mu"
    cli_args = [
        "--input", str(input_path),
        "--per-element-output", str(tmp_path / "per_element.tsv.gz"),
        "--per-guide-output", str(tmp_path / "per_guide.tsv.gz"),
        "--output-mudata", str(output_mudata),
    ]
    if test_all_pairs:
        cli_args.append("--test-all-pairs")
    args = adapter.build_parser().parse_args(cli_args)
    adapter.run_pipeline_adapter(args)

    assert output_mudata.exists()
    with h5py.File(output_mudata, "r") as f:
        x_after = f["mod/gene/X"][:].copy()
    assert np.array_equal(x_before, x_after), "X data should be byte-identical after a /uns-only patch"

    result = mu.read_h5mu(output_mudata)
    assert "per_element_results" in result.uns
    assert "per_guide_results" in result.uns
    assert f"{analysis_prefix}_per_element_results" in result.uns
    assert f"{analysis_prefix}_per_guide_results" in result.uns
    assert list(pd.DataFrame(result.uns["per_element_results"])["gene_id"]) == ["GENE1", "GENE2"]
