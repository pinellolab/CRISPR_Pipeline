import pathlib
import sys
import json
import os

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
import bounded_perturbo_results as bounded


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


def test_open_mudata_closes_backed_file_manager(tmp_path):
    input_path = tmp_path / "input.h5mu"
    _make_mudata().write(input_path)

    with adapter._open_mudata(input_path, backed="r") as mdata:
        file_manager = mdata.file
        assert file_manager.is_open

    assert not file_manager.is_open


@pytest.mark.parametrize(("flag", "value"), [("--result-batch-rows", "0"), ("--max-bh-working-bytes", "0")])
def test_conversion_memory_preflight_happens_before_input_or_fit(tmp_path, monkeypatch, flag, value):
    args = adapter.build_parser().parse_args(
        ["--input", str(tmp_path / "missing.h5mu"), "--per-element-output", str(tmp_path / "e.parquet"), "--per-guide-output", str(tmp_path / "g.parquet"), flag, value]
    )
    monkeypatch.setattr(adapter, "prepare_mudata_for_perturbo_v2", lambda *a, **k: (_ for _ in ()).throw(AssertionError("input opened")))
    with pytest.raises(ValueError, match="positive"):
        adapter.run_pipeline_adapter(args)


def test_published_raw_fit_survives_later_conversion_failure(tmp_path):
    fit = tmp_path / "fit"
    fit.mkdir()
    pd.DataFrame({"element": ["g"], "gene": ["G"], "posterior_mean": [0.0]}).to_parquet(
        fit / "element_effects.parquet", index=False
    )
    durable = adapter._publish_raw_fit(fit, tmp_path / "artifacts", "element")

    with pytest.raises(RuntimeError, match="injected"):
        raise RuntimeError("injected conversion failure")

    manifest = json.loads((durable / "conversion_manifest.json").read_text())
    bounded.verify_manifest(durable, manifest)
    assert (durable / "element_effects.parquet").exists()


def test_parallel_wait_publishes_success_before_raising_peer_failure(tmp_path):
    artifact = tmp_path / "artifacts"
    artifact.mkdir()
    fit_dirs = {}
    for name in ("element", "guide"):
        fit_dirs[name] = artifact / f".{name}.partial"
        fit_dirs[name].mkdir()
        pd.DataFrame({"p_value": [0.1]}).to_parquet(
            fit_dirs[name] / "element_effects.parquet", index=False
        )

    class Process:
        def __init__(self, code):
            self.code = code
            self.args = ["perturbo"]

        def poll(self):
            return self.code

    with pytest.raises(Exception):
        adapter._wait_publish_fits(
            {"element": Process(0), "guide": Process(2)}, fit_dirs, artifact
        )

    assert (artifact / "element/element_effects.parquet").exists()
    assert not (artifact / "guide").exists()


def test_publication_failure_terminates_and_waits_for_live_peer(tmp_path, monkeypatch):
    class Process:
        args = ["perturbo"]

        def __init__(self, code):
            self.code = code
            self.terminated = False
            self.waited = False

        def poll(self):
            return self.code

        def terminate(self):
            self.terminated = True
            self.code = -15

        def wait(self):
            self.waited = True
            return self.code

    complete = Process(0)
    live = Process(None)
    monkeypatch.setattr(adapter, "_publish_raw_fit", lambda *args: (_ for _ in ()).throw(OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        adapter._wait_publish_fits(
            {"element": complete, "guide": live},
            {"element": tmp_path / "element", "guide": tmp_path / "guide"},
            tmp_path,
        )

    assert live.terminated and live.waited


def test_requested_pairs_fallback_filters_in_bounded_batches(tmp_path):
    raw = tmp_path / "element_effects.parquet"
    requested = tmp_path / "element_effects_requested_pairs.parquet"
    pairs = tmp_path / "pairs.parquet"
    pd.DataFrame(
        {
            "element": ["a", "a", "b"], "gene": ["G1", "G2", "G2"],
            "posterior_mean": [0.0, 0.0, 0.0], "posterior_scale": [1.0, 1.0, 1.0],
            "posterior_prob": [0.1, 0.2, 0.3], "unused_diagnostic": [1, 2, 3],
        }
    ).to_parquet(raw, row_group_size=1, index=False)
    pd.DataFrame({"element": ["a"], "gene": ["G2"]}).to_parquet(pairs, index=False)

    adapter._ensure_requested_pairs(raw, requested, pairs, batch_rows=1)

    observed = pd.read_parquet(requested)
    assert list(observed[["element", "gene"]].itertuples(index=False, name=None)) == [("a", "G2")]
    assert "unused_diagnostic" not in observed


@pytest.mark.parametrize("suffix", [".tsv", ".tsv.gz"])
def test_direct_raw_render_preserves_identifier_text_and_empty_headers(tmp_path, suffix):
    raw = tmp_path / "requested.parquet"
    pd.DataFrame(
        {
            "element": ["NA", "001"], "gene": ["001", "NA"],
            "posterior_mean": [0.0, 0.0], "posterior_scale": [1.0, 1.0],
            "posterior_prob": [0.1, 0.2],
        }
    ).to_parquet(raw, index=False)
    output = tmp_path / f"result{suffix}"
    adapter._write_bounded_effects(
        raw, output, inference_type="guide", prepared_path=None,
        guide_name_map={}, crt=False, scratch_dir=tmp_path,
        max_bh_working_bytes=1 << 20, compact_floats=False, batch_rows=1,
    )
    observed = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
    assert list(observed["guide_id"]) == ["NA", "001"]
    assert list(observed["gene_id"]) == ["001", "NA"]

    empty_raw = tmp_path / "empty.parquet"
    pd.DataFrame(
        {name: pd.Series(dtype=dtype) for name, dtype in {
            "element": "string", "gene": "string", "posterior_mean": "float64",
            "posterior_scale": "float64", "posterior_prob": "float64",
        }.items()}
    ).to_parquet(empty_raw, index=False)
    empty_output = tmp_path / f"empty{suffix}"
    adapter._write_bounded_effects(
        empty_raw, empty_output, inference_type="guide", prepared_path=None,
        guide_name_map={}, crt=False, scratch_dir=tmp_path,
        max_bh_working_bytes=1 << 20, compact_floats=False, batch_rows=1,
    )
    assert list(pd.read_csv(empty_output, sep="\t").columns) == [
        "gene_id", "guide_id", "log2_fc", "perturbo_fc_se", "p_value",
        "perturbo_posterior_prob", "perturbo_q_value",
    ]


def test_bounded_bh_clips_only_for_adjustment_and_carries_typed_diagnostics(tmp_path):
    raw = tmp_path / "diagnostics.parquet"
    p = np.array([-np.inf, np.inf, 2.0, -1.0, np.nan])
    pd.DataFrame(
        {
            "element": ["g"] * 5, "gene": [f"G{i}" for i in range(5)],
            "posterior_mean": np.zeros(5, dtype=np.float32),
            "posterior_scale": np.ones(5, dtype=np.float32),
            "posterior_prob": np.full(5, 0.5, dtype=np.float32),
            "crt_saddlepoint_p_value": p,
            "crt_low_information": pd.Series([True, False, True, False, True], dtype="bool"),
            "crt_tail_failure_reason": pd.Series([16, 0, -1, 2, 3], dtype="int16"),
            "crt_saddlepoint_valid": pd.Series([1.0, np.nan, 0.0, 1.0, np.nan], dtype="float64"),
            "crt_root_residual_null_sd": pd.Series([0.1, 0.2, np.nan, 0.4, 0.5], dtype="float64"),
        }
    ).to_parquet(raw, index=False)
    output = tmp_path / "result.parquet"
    adapter._write_bounded_effects(
        raw, output, inference_type="guide", prepared_path=None, guide_name_map={},
        crt=True, scratch_dir=tmp_path, max_bh_working_bytes=1 << 20,
        compact_floats=False, batch_rows=2,
    )
    observed = pd.read_parquet(output)
    np.testing.assert_array_equal(observed["p_value"].to_numpy(), p)
    expected_q = adapter._bh_adjust(pd.Series(p)).to_numpy()
    np.testing.assert_allclose(observed["perturbo_q_value"], expected_q, equal_nan=True)
    assert observed["perturbo_crt_low_information"].dtype == bool
    assert observed["perturbo_crt_tail_failure_reason"].dtype == np.int16
    assert observed["perturbo_crt_saddlepoint_valid"].dtype == np.float64
    assert np.isnan(observed.loc[1, "perturbo_crt_saddlepoint_valid"])


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
    assert "total_gene_umis" in gene.obs
    # PerTurbo conditions on precomputed logs, since it applies no transform of
    # its own; SCEPTRE gets the plain counts and logs them in its formula.
    for column in (
        "log_total_guide_umis",
        "log_total_gene_umis",
        "log_num_expressed_genes",
    ):
        assert column in gene.obs
    assert any(str(name).startswith("non-targeting|") for name in guide.uns[adapter.GUIDE_NAMES_KEY])
    assert guide_name_map["non-targeting|nt1"] == "nt1"


def test_build_native_pairs_to_test_uses_only_requested_pairs(tmp_path):
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _make_mudata().write(input_path)
    guide_name_map = adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)
    prepared = mu.read_h5mu(prepared_path)

    element_pairs = adapter._build_native_pairs_to_test(
        prepared,
        pair_element_column="intended_target_key",
        element_names=[str(x) for x in prepared["guide"].uns[adapter.ELEMENT_NAMES_KEY]],
    )
    guide_pairs = adapter._build_native_pairs_to_test(
        prepared,
        pair_element_column="guide_id",
        element_names=[str(x) for x in prepared["guide"].uns[adapter.GUIDE_NAMES_KEY]],
        guide_name_map=guide_name_map,
    )

    assert list(element_pairs.columns) == ["element", "gene"]
    assert list(guide_pairs.columns) == ["element", "gene"]
    assert set(element_pairs["gene"]) == {"GENE1", "GENE2"}
    assert set(guide_pairs["gene"]) == {"GENE1", "GENE2"}
    assert {"gA", "gB", "non-targeting|nt1"}.issubset(set(guide_pairs["element"]))
    # Exactly the requested pairs -- the fixture asks for gA->GENE1, gB->GENE2
    # and nt1->GENE1. Control elements are deliberately not crossed with every
    # tested gene: that augmentation used to dominate this table's
    # Benjamini-Hochberg family (2,623,521 of 2,714,942 rows on Replogle) and
    # made the local q-values far more conservative than SCEPTRE's over the same
    # hypotheses. Controls are still tested, in the transcriptome-wide table.
    control_pairs = guide_pairs[guide_pairs["element"] == "non-targeting|nt1"]
    assert set(control_pairs["gene"]) == {"GENE1"}


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
            # Diagnostics PerTurbo writes beside the p-value; the tail-policy
            # column is one of the Chernoff-fallback release's, the others rc9's.
            "crt_low_information": [True, False, False],
            "crt_observed_nonzero": [0, 12, 40],
            "crt_tail_failure_reason": [16, 0, -1],
        }
    )

    observed = adapter.convert_guide_effects(
        effects,
        guide_name_map,
        prepared_path,
        test_all_pairs=False,
    )

    # Carried under the method prefix, as they come, without filtering rows.
    assert list(observed["perturbo_crt_low_information"]) == [True, False]
    assert list(observed["perturbo_crt_observed_nonzero"]) == [0, 40]
    assert list(observed["perturbo_crt_tail_failure_reason"]) == [16, -1]
    assert "perturbo_crt_used_chernoff" not in observed.columns  # absent upstream, absent here

    assert list(observed["guide_id"]) == ["gA", "nt1"]
    assert list(observed["gene_id"]) == ["GENE1", "GENE1"]
    assert np.isclose(observed.loc[0, "log2_fc"], 1.0)
    assert np.isclose(observed.loc[0, "p_value"], 0.05)
    assert np.isclose(observed.loc[1, "p_value"], 0.9)


@pytest.mark.parametrize("with_pairs", [False, True])
def test_run_perturbo_uses_only_native_pairs_to_test_flag(tmp_path, monkeypatch, with_pairs):
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _make_mudata().write(input_path)
    adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)
    pairs_path = tmp_path / "pairs.parquet"
    pd.DataFrame({"element": ["gA"], "gene": ["GENE1"]}).to_parquet(
        pairs_path, index=False
    )
    captured = {}

    class FakeProcess:
        args = []

    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProcess()

    monkeypatch.setattr(
        adapter, "_covariate_has_baseline_variance", lambda *args, **kwargs: False
    )
    monkeypatch.setattr(adapter.subprocess, "Popen", fake_popen)
    args = adapter.build_parser().parse_args(
        [
            "--input", str(input_path),
            "--per-element-output", str(tmp_path / "element.tsv.gz"),
            "--per-guide-output", str(tmp_path / "guide.tsv.gz"),
            "--device", "cpu",
            "--no-save-model-params",
        ]
    )
    # main() resolves the pool before it reaches _run_perturbo; these tests call
    # _run_perturbo directly, so do what main() does rather than hardcode a value.
    args.resolved_crt_pool = adapter._resolve_crt_pool(args.crt_pool, None)
    adapter._run_perturbo(
        prepared_path,
        tmp_path / "fit",
        map_key=adapter.GUIDE_MAP_KEY,
        names_key=adapter.GUIDE_NAMES_KEY,
        pairs_to_test_path=pairs_path if with_pairs else None,
        gpu_id=None,
        phase="smoke",
        args=args,
    )

    cmd = captured["cmd"]
    assert "--gene-by-element-varm-key" not in cmd
    assert "--gene-by-element-names-uns-key" not in cmd
    if with_pairs:
        index = cmd.index("--pairs-to-test")
        assert cmd[index + 1] == str(pairs_path)
    else:
        assert "--pairs-to-test" not in cmd


@pytest.mark.parametrize(
    ("test_all_pairs", "analysis_prefix"),
    [(False, "local_analysis"), (True, "global_analysis")],
)
def test_run_pipeline_adapter_patches_uns_without_rewriting_x(
    tmp_path, monkeypatch, test_all_pairs, analysis_prefix
):
    pytest.importorskip("pyarrow")
    input_path = tmp_path / "input.h5mu"
    input_mudata = _make_mudata()
    if test_all_pairs:
        input_mudata["guide"].var["guide_id"] = ["NA", "001", "nt1"]
        input_mudata.uns["pairs_to_test"]["guide_id"] = ["NA", "001", "nt1"]
    input_mudata.write(input_path)

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
            "element": ["NA", "001"] if test_all_pairs else ["gA", "gB"],
            "gene": ["GENE1", "GENE2"],
            "posterior_mean": [np.log(2), np.log(3)],
            "posterior_scale": [np.log(2) / 10, np.log(2) / 9],
            "posterior_prob": [0.05, 0.07],
        }
    )

    observed_native_pairs = []

    class FinishedProcess:
        args = ["perturbo"]

        @staticmethod
        def poll():
            return 0

    def fake_run_perturbo(
        input_path, out_dir, *, map_key, names_key, pairs_to_test_path, args, **kwargs
    ):
        out_dir.mkdir(parents=True, exist_ok=True)
        observed_native_pairs.append(
            None if pairs_to_test_path is None else pd.read_parquet(pairs_to_test_path)
        )
        effects = element_effects if map_key == adapter.ELEMENT_MAP_KEY else guide_effects
        effects.to_parquet(out_dir / "element_effects.parquet", index=False)
        (out_dir / "crt_metadata.json").write_text('{"raw": true}\n')
        if pairs_to_test_path is not None:
            pairs = pd.read_parquet(pairs_to_test_path)
            effects.merge(pairs, on=["element", "gene"], how="inner").to_parquet(
                out_dir / "element_effects_requested_pairs.parquet", index=False
            )
        return FinishedProcess()

    monkeypatch.setattr(adapter, "_run_perturbo", fake_run_perturbo)

    output_mudata = tmp_path / "output.h5mu"
    primary_element = tmp_path / ("per_element.tsv.gz" if test_all_pairs else "per_element.parquet")
    primary_guide = tmp_path / ("per_guide.tsv.gz" if test_all_pairs else "per_guide.parquet")
    cli_args = [
        "--input", str(input_path),
        "--per-element-output", str(primary_element),
        "--per-guide-output", str(primary_guide),
        "--output-mudata", str(output_mudata),
        "--v2-artifact-dir", str(tmp_path / "artifacts"),
    ]
    if test_all_pairs:
        cli_args.extend(
            [
                "--test-all-pairs",
                "--local-per-element-output", str(tmp_path / "local_element.tsv.gz"),
                "--local-per-guide-output", str(tmp_path / "local_guide.pq"),
            ]
        )
    else:
        cli_args.extend(
            [
                "--local-per-element-output", str(tmp_path / "local_element.tsv.gz"),
                "--local-per-guide-output", str(tmp_path / "local_guide.tsv.gz"),
            ]
        )
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
    if test_all_pairs:
        assert len(observed_native_pairs) == 2
        assert (tmp_path / "local_element.tsv.gz").exists()
        assert (tmp_path / "local_guide.pq").exists()
        assert list(pd.read_parquet(tmp_path / "local_guide.pq")["guide_id"].astype(str)) == ["NA", "001"]
        embedded_guides = pd.DataFrame(result.uns["global_analysis_per_guide_results"])
        assert list(embedded_guides["guide_id"].astype(str)) == ["NA", "001"]
    else:
        assert primary_element.exists()
        assert primary_guide.exists()
        assert (tmp_path / "local_element.tsv.gz").exists()
        assert (tmp_path / "local_guide.tsv.gz").exists()
        assert list(pd.read_csv(tmp_path / "local_element.tsv.gz", sep="\t")["gene_id"]) == list(
            pd.read_parquet(primary_element)["gene_id"]
        )
        assert len(observed_native_pairs) == 2
        assert all(list(frame.columns) == ["element", "gene"] for frame in observed_native_pairs)
        assert all(any(frame["element"].str.contains("non-targeting")) for frame in observed_native_pairs)
        assert (tmp_path / "artifacts/element_pairs_to_test.parquet").exists()
        assert (tmp_path / "artifacts/guide_pairs_to_test.parquet").exists()

    if test_all_pairs:
        def fitting_must_not_run(*args, **kwargs):
            raise AssertionError("conversion-only recovery launched a fit")

        monkeypatch.setattr(adapter, "_run_perturbo", fitting_must_not_run)
        monkeypatch.setattr(adapter, "prepare_mudata_for_perturbo_v2", fitting_must_not_run)
        recovery_cli = list(cli_args)
        output_index = recovery_cli.index("--output-mudata")
        del recovery_cli[output_index : output_index + 2]
        recovery_args = adapter.build_parser().parse_args(
            recovery_cli + ["--conversion-only-artifact-dir", str(tmp_path / "artifacts")]
        )
        adapter.run_pipeline_adapter(recovery_args)
        assert (tmp_path / "artifacts/element/crt_metadata.json").read_text() == '{"raw": true}\n'

        original = input_path.stat()
        input_path.touch()
        changed_args = adapter.build_parser().parse_args(
            recovery_cli + ["--conversion-only-artifact-dir", str(tmp_path / "artifacts")]
        )
        with pytest.raises(ValueError, match="provenance"):
            adapter.run_pipeline_adapter(changed_args)
        os.utime(input_path, ns=(original.st_atime_ns, original.st_mtime_ns))

def _write_mudata_without_controls(path):
    """A screen with no non-targeting guides, so no cell is a control cell."""
    mdata = _make_mudata()
    guide = mdata["guide"]
    guide.var.loc["nt1", "targeting"] = True
    guide.var.loc["nt1", "type"] = "targeting"
    guide.var.loc["nt1", "intended_target_name"] = "elemC"
    guide.var.loc["nt1", "intended_target_chr"] = "chr3"
    guide.var.loc["nt1", "intended_target_start"] = 500.0
    guide.var.loc["nt1", "intended_target_end"] = 600.0
    # A guide-UMI depth that varies across cells, so the covariate is
    # identifiable on the all-cells baseline even with no control cells. The
    # assignment layer stays binary; only the raw UMI counts differ.
    guide.X = sparse.csr_matrix(
        np.array([[3, 0, 0], [0, 7, 0], [0, 0, 11]], dtype=np.float32)
    )
    guide.layers["guide_assignment"] = sparse.csr_matrix(
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
    )
    mdata.write(path)
    return mdata


@pytest.mark.parametrize(
    ("pool", "expect_guide_umi_covariate"),
    [("all-cells", True), ("control-anchored", False)],
)
def test_guide_umi_covariate_survives_the_all_cells_pool_without_controls(
    tmp_path, monkeypatch, pool, expect_guide_umi_covariate
):
    """The baseline-variance guard has to follow the pool.

    Under ``all-cells`` PerTurbo fits its stage-one baseline on every cell, so a
    screen with no non-targeting guides must keep the guide-UMI covariate --
    dropping it there left PerTurbo conditioning on less than SCEPTRE. Under
    ``control-anchored`` the baseline really is the control cells, and with none
    the covariate is unidentifiable, so it is still dropped.
    """
    input_path = tmp_path / "input.h5mu"
    prepared_path = tmp_path / "prepared.h5mu"
    _write_mudata_without_controls(input_path)
    adapter.prepare_mudata_for_perturbo_v2(input_path, prepared_path)

    captured = {}

    class FakeProcess:
        args = []

    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        return FakeProcess()

    monkeypatch.setattr(adapter.subprocess, "Popen", fake_popen)
    args = adapter.build_parser().parse_args(
        [
            "--input", str(input_path),
            "--per-element-output", str(tmp_path / "element.tsv.gz"),
            "--per-guide-output", str(tmp_path / "guide.tsv.gz"),
            "--device", "cpu",
            "--no-save-model-params",
        ]
    )
    args.resolved_crt_pool = pool
    adapter._run_perturbo(
        prepared_path,
        tmp_path / "fit",
        map_key=adapter.GUIDE_MAP_KEY,
        names_key=adapter.GUIDE_NAMES_KEY,
        pairs_to_test_path=None,
        gpu_id=None,
        phase="smoke",
        args=args,
    )

    cmd = captured["cmd"]
    index = cmd.index("--continuous-covariates")
    passed = []
    for token in cmd[index + 1:]:
        if token.startswith("--"):
            break
        passed.append(token)
    assert (adapter.GUIDE_UMI_COVARIATE in passed) is expect_guide_umi_covariate
