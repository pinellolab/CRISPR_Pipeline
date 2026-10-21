"""The calibration and AUPRC report over matched non-targeting pairs.

Inputs are synthesised from the column contracts of the prepare step, the
SCEPTRE chunk merge and the PerTurbo adapter, with uniform p-values so the
rates land near nominal.
"""

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import ntc_calibration_report as report

GENES = [f"G{i}" for i in range(1, 7)]
BUCKETS = {1: ["nt1", "nt2"], 2: ["nt3", "nt4"]}
TARGETING = {"g1": "G1", "g2": "G1", "g3": "G2", "g4": "G2"}


def _inputs(tmp_path, *, with_main_sceptre=True, seed=11):
    rng = np.random.default_rng(seed)
    sets = {}
    for bucket in BUCKETS:
        for gene in GENES[:4]:
            sets.setdefault((bucket, gene), set()).add("global")
        for gene in ("G1", "G2"):
            sets.setdefault((bucket, gene), set()).add("direct_target")
    for bucket, gene in ((1, "G1"), (1, "G3"), (2, "G2"), (2, "G5")):
        sets.setdefault((bucket, gene), set()).add("local")
    rows = []
    for (bucket, gene), labels in sorted(sets.items()):
        for guide in BUCKETS[bucket]:
            rows.append({
                "guide_id": guide, "gene_id": gene, "intended_target_name": f"ntc_calibration|{bucket}",
                "intended_target_chr": "", "intended_target_start": pd.NA, "intended_target_end": pd.NA,
                "intended_target_key": f"name_only::ntc_calibration|{bucket}", "pair_type": "ntc_calibration",
                "pair_set": "+".join(sorted(labels)), "ntc_bucket": bucket, "substituted_element": None,
            })
    pairs = pd.DataFrame(rows)
    paths = {"pairs": tmp_path / "pairs.tsv.gz"}
    pairs.to_csv(paths["pairs"], sep="\t", index=False, compression="gzip")

    guide_targets = pd.DataFrame(
        [{"guide_id": g, "intended_target_name": t, "intended_target_key": f"coords::{t}::chr1::1::2", "targeting": True, "role": "targeting", "ntc_bucket": pd.NA} for g, t in TARGETING.items()]
        + [{"guide_id": g, "intended_target_name": f"non-targeting|{b}", "intended_target_key": f"coords::non-targeting|{b}::non-targeting::{b}::{b}", "targeting": False, "role": "ntc", "ntc_bucket": b} for b, gs in BUCKETS.items() for g in gs]
    )
    paths["guide_targets"] = tmp_path / "guide_targets.tsv"
    guide_targets.to_csv(paths["guide_targets"], sep="\t", index=False)
    direct = pd.DataFrame([{"guide_id": g, "gene_id": t} for g, t in TARGETING.items()])
    paths["direct_targets"] = tmp_path / "direct.tsv"
    direct.to_csv(paths["direct_targets"], sep="\t", index=False)

    def sceptre_guide(keys):
        return pd.DataFrame({"gene_id": [k[1] for k in keys], "guide_id": [k[0] for k in keys],
                             "p_value": rng.uniform(size=len(keys)), "q_value": np.nan,
                             "fold_change": 1.0, "se_fold_change": 0.1, "log2_fc": rng.normal(size=len(keys))})

    def sceptre_element(keys):
        return pd.DataFrame({"gene_id": [k[1] for k in keys], "intended_target_name": [k[0] for k in keys],
                             "intended_target_chr": "", "intended_target_start": pd.NA, "intended_target_end": pd.NA,
                             "fold_change": 1.0, "se_fold_change": 0.1, "log2_fc": rng.normal(size=len(keys)),
                             "p_value": rng.uniform(size=len(keys)), "q_value": np.nan})

    pass_guide = sceptre_guide(list(zip(pairs["guide_id"], pairs["gene_id"])))
    pass_element = sceptre_element([(f"ntc_calibration|{b}", g) for (b, g) in sets])
    paths["sceptre_ntc_per_guide"] = tmp_path / "ntc_sceptre_per_guide_output.tsv.gz"
    paths["sceptre_ntc_per_element"] = tmp_path / "ntc_sceptre_per_element_output.tsv.gz"
    pass_guide.to_csv(paths["sceptre_ntc_per_guide"], sep="\t", index=False, compression="gzip")
    pass_element.to_csv(paths["sceptre_ntc_per_element"], sep="\t", index=False, compression="gzip")

    all_guides = list(TARGETING) + [g for gs in BUCKETS.values() for g in gs]
    global_guide = pd.DataFrame([(g, gene) for g in all_guides for gene in GENES], columns=["guide_id", "gene_id"])
    global_guide["log2_fc"] = rng.normal(size=len(global_guide))
    global_guide["perturbo_fc_se"] = 0.1
    global_guide["p_value"] = rng.uniform(size=len(global_guide))
    # Direct targets are strong hits for PerTurbo.
    direct_mask = global_guide.apply(lambda r: TARGETING.get(r["guide_id"]) == r["gene_id"], axis=1)
    global_guide.loc[direct_mask, "p_value"] = 1e-6
    global_guide["perturbo_q_value"] = np.where(global_guide["guide_id"].isin(TARGETING), 0.5, np.nan)
    global_guide["perturbo_crt_low_information"] = rng.random(len(global_guide)) < 0.2
    elements = ["G1", "G2", "non-targeting|1", "non-targeting|2"]
    global_element = pd.DataFrame([(e, gene) for e in elements for gene in GENES], columns=["intended_target_name", "gene_id"])
    global_element["intended_target_chr"] = ""
    global_element["intended_target_start"] = pd.NA
    global_element["intended_target_end"] = pd.NA
    global_element["log2_fc"] = rng.normal(size=len(global_element))
    global_element["perturbo_fc_se"] = 0.1
    global_element["p_value"] = rng.uniform(size=len(global_element))
    global_element.loc[global_element["intended_target_name"] == global_element["gene_id"], "p_value"] = 1e-6
    global_element["perturbo_q_value"] = np.nan
    global_element["perturbo_crt_low_information"] = rng.random(len(global_element)) < 0.2
    paths["perturbo_global_per_guide"] = tmp_path / "perturbo_global_analysis_per_guide_output.parquet"
    paths["perturbo_global_per_element"] = tmp_path / "perturbo_global_analysis_per_element_output.parquet"
    global_guide.to_parquet(paths["perturbo_global_per_guide"], index=False)
    global_element.to_parquet(paths["perturbo_global_per_element"], index=False)

    if with_main_sceptre:
        main_guide = sceptre_guide([(g, gene) for g, t in TARGETING.items() for gene in (t, "G3")])
        main_guide.loc[main_guide.apply(lambda r: TARGETING[r["guide_id"]] == r["gene_id"], axis=1), "p_value"] = 1e-5
        main_element = sceptre_element([(t, gene) for t in ("G1", "G2") for gene in (t, "G3")])
        main_element.loc[main_element["intended_target_name"] == main_element["gene_id"], "p_value"] = 1e-5
        paths["sceptre_local_per_guide"] = tmp_path / "main_sceptre_per_guide_output.tsv.gz"
        paths["sceptre_local_per_element"] = tmp_path / "main_sceptre_per_element_output.tsv.gz"
        main_guide.to_csv(paths["sceptre_local_per_guide"], sep="\t", index=False, compression="gzip")
        main_element.to_csv(paths["sceptre_local_per_element"], sep="\t", index=False, compression="gzip")

    design = {"status": "ok", "regime": "complement", "analysis_moi": "high", "reference_pool_size": 100,
              "self_inclusion_ratio_by_bucket": {"1": 0.05, "2": 0.04}}
    paths["summary_in"] = tmp_path / "design.json"
    paths["summary_in"].write_text(json.dumps(design))
    paths["inference_mode"] = tmp_path / "inference_mode.json"
    paths["inference_mode"].write_text(json.dumps({"analysis_moi": "high", "sceptre_control_group": "complement", "perturbo_crt_pool": "all-cells"}))
    per_guide = pd.DataFrame({"guide_id": all_guides, "role": ["targeting"] * 4 + ["ntc"] * 4,
                              "n_cells": [30, 28, 35, 31, 29, 33, 27, 30], "n_cells_singleton": [28, 27, 33, 30, 29, 33, 27, 30]})
    per_element = pd.DataFrame({"intended_target_key": ["k1", "k2", "k3", "k4"], "element_name": ["G1", "G2", "non-targeting|1", "non-targeting|2"],
                                "role": ["targeting", "targeting", "ntc", "ntc"], "n_guides": 2, "ntc_bucket": [pd.NA, pd.NA, 1, 2], "n_cells": [58, 66, 62, 57]})
    paths["cell_counts_per_guide"] = tmp_path / "cells_guide.tsv"
    paths["cell_counts_per_element"] = tmp_path / "cells_element.tsv"
    per_guide.to_csv(paths["cell_counts_per_guide"], sep="\t", index=False)
    per_element.to_csv(paths["cell_counts_per_element"], sep="\t", index=False)
    return paths


def _argv(paths, outdir, skipped=False):
    argv = []
    if skipped:
        argv.append("--skipped")
        wanted = ("summary_in", "cell_counts_per_guide", "cell_counts_per_element")
    else:
        wanted = paths.keys()
    for name in wanted:
        argv += [f"--{name.replace('_', '-')}", str(paths[name])]
    return argv + ["--outdir", str(outdir)]


def test_report_writes_summaries_plots_and_joined_tables(tmp_path):
    paths = _inputs(tmp_path)
    outdir = tmp_path / "calibration"
    assert report.main(_argv(paths, outdir)) == 0

    for name in (
        "ntc_calibration_per_guide.tsv.gz", "ntc_calibration_per_element.tsv.gz", "ntc_calibration_summary.tsv",
        "ntc_calibration_summary.json", "ntc_method_agreement.tsv", "ntc_aucpr_summary.tsv", "ntc_cell_count_summary.tsv",
        "ntc_calibration_notes.txt", "ntc_calibration_qq_guide.png", "ntc_calibration_qq_element.png",
        "ntc_calibration_pvalue_hist_guide.png", "ntc_method_scatter_guide.png", "ntc_perturbo_global_qq.png",
        "ntc_aucpr_curves_guide.png", "ntc_aucpr_curves_element.png", "ntc_cells_per_guide.png", "ntc_cells_per_element.png",
    ):
        assert (outdir / name).exists(), name

    summary = pd.read_csv(outdir / "ntc_calibration_summary.tsv", sep="\t")
    expected = {(lvl, ps, m) for lvl in ("guide", "element") for ps in ("local", "global", "direct_target") for m in ("sceptre", "perturbo")}
    expected |= {(lvl, "perturbo_global_all", "perturbo") for lvl in ("guide", "element")}
    assert set(map(tuple, summary[["level", "pair_set", "method"]].to_numpy())) == expected
    assert summary["n_tested"].gt(0).all()
    assert summary["frac_p_lt_0.05"].between(0, 1).all()
    assert summary.loc[summary["method"] == "perturbo", "n_low_information"].notna().all()

    joined = pd.read_csv(outdir / "ntc_calibration_per_guide.tsv.gz", sep="\t")
    assert joined[["sceptre_p_value", "perturbo_p_value"]].notna().all().all()
    assert "n_cells_guide" in joined.columns and "perturbo_crt_low_information" in joined.columns
    assert not any(c.endswith("q_value") for c in joined.columns)

    agreement = pd.read_csv(outdir / "ntc_method_agreement.tsv", sep="\t")
    assert agreement["n_both_tested"].gt(0).all()
    assert agreement["pearson_log10p"].notna().any()

    auprc = pd.read_csv(outdir / "ntc_aucpr_summary.tsv", sep="\t")
    assert set(map(tuple, auprc[["level", "method"]].to_numpy())) == {(l, m) for l in ("guide", "element") for m in ("sceptre", "perturbo")}
    assert auprc["auprc"].between(0, 1).all() and auprc["auroc"].between(0, 1).all()
    # Direct targets were planted as strong hits, so both methods separate them.
    assert auprc["auroc"].gt(0.9).all()

    notes = (outdir / "ntc_calibration_notes.txt").read_text()
    assert "Self-inclusion" in notes and "5.0%" in notes


def test_report_tolerates_a_missing_main_sceptre_table(tmp_path):
    paths = _inputs(tmp_path, with_main_sceptre=False)
    outdir = tmp_path / "calibration"
    assert report.main(_argv(paths, outdir)) == 0
    auprc = pd.read_csv(outdir / "ntc_aucpr_summary.tsv", sep="\t")
    assert set(auprc["method"]) == {"perturbo"}
    assert "AUPRC is not reported" in (outdir / "ntc_calibration_notes.txt").read_text()


def test_untested_pairs_count_as_untested_not_as_hits(tmp_path):
    paths = _inputs(tmp_path)
    frame = pd.read_csv(paths["sceptre_ntc_per_guide"], sep="\t")
    frame.loc[:3, "p_value"] = np.nan
    frame.to_csv(paths["sceptre_ntc_per_guide"], sep="\t", index=False, compression="gzip")
    outdir = tmp_path / "calibration"
    assert report.main(_argv(paths, outdir)) == 0
    summary = pd.read_csv(outdir / "ntc_calibration_summary.tsv", sep="\t")
    sceptre_guide = summary[(summary["level"] == "guide") & (summary["method"] == "sceptre")]
    assert (sceptre_guide["n_tested"] <= sceptre_guide["n_pairs"]).all()
    assert (sceptre_guide["n_tested"] < sceptre_guide["n_pairs"]).any()


def test_skip_mode_writes_the_note_and_the_cell_counts_only(tmp_path):
    paths = _inputs(tmp_path)
    paths["summary_in"].write_text(json.dumps({"status": "skipped", "reason": "the screen has no non-targeting guides"}))
    outdir = tmp_path / "calibration"
    assert report.main(_argv(paths, outdir, skipped=True)) == 0
    assert "no non-targeting guides" in (outdir / "ntc_calibration_skipped.txt").read_text()
    assert (outdir / "ntc_cell_count_summary.tsv").exists()
    assert (outdir / "ntc_cells_per_guide.png").exists()
    assert not (outdir / "ntc_calibration_summary.tsv").exists()


def test_cell_count_summary_warns_when_ntc_cells_are_off_by_more_than_twofold():
    per_guide = pd.DataFrame({"guide_id": ["a", "b", "c", "d"], "role": ["targeting", "targeting", "ntc", "ntc"], "n_cells": [100, 120, 10, 12]})
    per_element = pd.DataFrame({"role": ["targeting", "ntc"], "n_cells": [200, 190]})
    summary, warnings = report.cell_count_summary(per_guide, per_element)
    assert len(warnings) == 1 and warnings[0].startswith("guide level")
    assert summary.loc[summary["level"] == "guide", "warning"].str.len().gt(0).all()
    assert summary.loc[summary["level"] == "element", "warning"].eq("").all()
