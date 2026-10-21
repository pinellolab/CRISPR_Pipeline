#!/usr/bin/env python
"""Calibration and AUPRC of SCEPTRE and PerTurbo on matched non-targeting pairs.

Joins the SCEPTRE calibration pass to PerTurbo's global tables on the pair
table the prepare step wrote, and reports, per level (guide, element), per pair
set and per method: how the non-targeting p-values are distributed against the
uniform (QQ plots, empirical rates at 0.01, 0.05 and 0.10, a KS distance), how
the two methods agree on the same pairs, and the direct-target-versus-NTC
AUPRC and AUROC the controls evaluation has always reported for PerTurbo, now
for both methods. The sanity check compares cells per guide and per element
between the non-targeting and targeting sets.

Only p-values are read. The q-values the reused processes wrote are over the
calibration family and are not reported.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from inference_target_matching import normalize_target_identifier  # noqa: E402
from result_table_io import read_result_table  # noqa: E402

CALIBRATION_PREFIX = "ntc_calibration|"
CONTROL_PREFIX = "non-targeting|"
P_FLOOR = 1e-12
THRESHOLDS = (0.01, 0.05, 0.10)
AUPRC_SEED = 42
LEVELS = ("guide", "element")
METHODS = ("sceptre", "perturbo")


# --- loading -----------------------------------------------------------------


def _bucket(name) -> int | None:
    text = str(name)
    if "|" not in text:
        return None
    try:
        return int(text.rsplit("|", 1)[1])
    except ValueError:
        return None


def _p_column(frame: pd.DataFrame, method: str) -> str | None:
    for name in (f"{method}_p_value", "p_value"):
        if name in frame.columns:
            return name
    return None


def _fc_column(frame: pd.DataFrame, method: str) -> str | None:
    for name in (f"{method}_log2_fc", "log2_fc"):
        if name in frame.columns:
            return name
    return None


def load_method_table(path, method: str, level: str) -> pd.DataFrame:
    """A result table reduced to its join keys, p-value, fold change and diagnostics."""
    frame = read_result_table(path)
    key = "guide_id" if level == "guide" else "intended_target_name"
    p_col, fc_col = _p_column(frame, method), _fc_column(frame, method)
    if p_col is None or key not in frame.columns or "gene_id" not in frame.columns:
        raise ValueError(f"{path}: expected {key}, gene_id and a p-value column for {method}")
    out = pd.DataFrame(
        {
            key: frame[key].astype(str).to_numpy(),
            "gene_id": frame["gene_id"].astype(str).to_numpy(),
            f"{method}_p_value": pd.to_numeric(frame[p_col], errors="coerce").to_numpy(dtype=float),
            f"{method}_log2_fc": pd.to_numeric(frame[fc_col], errors="coerce").to_numpy(dtype=float) if fc_col else np.nan,
        }
    )
    for column in frame.columns:
        if column.startswith("perturbo_crt_") and method == "perturbo":
            out[column] = frame[column].to_numpy()
    return out.drop_duplicates(subset=[key, "gene_id"], keep="first")


# --- statistics ----------------------------------------------------------------


def calibration_rows(p: np.ndarray, *, level: str, pair_set: str, method: str, low_information=None) -> dict:
    finite = np.isfinite(p)
    tested = p[finite]
    row = {
        "level": level, "pair_set": pair_set, "method": method,
        "n_pairs": int(p.size), "n_tested": int(tested.size),
    }
    for threshold in THRESHOLDS:
        row[f"frac_p_lt_{threshold:.2f}"] = float((tested < threshold).mean()) if tested.size else np.nan
    row["ks_statistic"] = float(stats.kstest(tested, "uniform").statistic) if tested.size else np.nan
    row["median_p"] = float(np.median(tested)) if tested.size else np.nan
    if low_information is not None:
        row["n_low_information"] = int(np.nansum(np.asarray(low_information, dtype=float) > 0))
    return row


def agreement_row(joined: pd.DataFrame, *, level: str, pair_set: str) -> dict:
    both = joined.dropna(subset=["sceptre_p_value", "perturbo_p_value"])
    row = {"level": level, "pair_set": pair_set, "n_both_tested": int(len(both))}
    if len(both) >= 3:
        s = -np.log10(np.clip(both["sceptre_p_value"], P_FLOOR, 1.0))
        t = -np.log10(np.clip(both["perturbo_p_value"], P_FLOOR, 1.0))
        row["pearson_log10p"] = float(np.corrcoef(s, t)[0, 1]) if s.std() > 0 and t.std() > 0 else np.nan
        row["spearman_p"] = float(stats.spearmanr(both["sceptre_p_value"], both["perturbo_p_value"]).statistic)
        a, b = both["sceptre_p_value"] < 0.05, both["perturbo_p_value"] < 0.05
        union = int((a | b).sum())
        row["jaccard_0.05"] = float((a & b).sum() / union) if union else np.nan
    else:
        row.update({"pearson_log10p": np.nan, "spearman_p": np.nan, "jaccard_0.05": np.nan})
    return row


def auprc_rows(positives: pd.DataFrame, negatives: pd.DataFrame, *, level: str, method: str, p_col: str) -> list[dict]:
    """The controls evaluation's metric, twice: direct-target rows against NTC
    rows on the same genes, scored by 1 - p, with the negatives drawn 1:1 as
    the evaluation has always done (``matched_1to1``) and with every available
    negative (``all``), which moves the precision-recall baseline to the real
    positive rate but removes the draw's sampling noise."""
    from sklearn.metrics import auc, precision_recall_curve, roc_curve

    pos = positives.dropna(subset=[p_col])
    neg = negatives.dropna(subset=[p_col])
    neg = neg[neg["gene_id"].isin(set(pos["gene_id"]))] if len(pos) else neg
    base = {"level": level, "method": method, "n_direct_target": int(len(pos)), "n_ntc_available": int(len(neg))}
    rows = []
    for variant in ("matched_1to1", "all"):
        row = {**base, "negatives": variant}
        if pos.empty or neg.empty:
            row.update({"n_ntc_used": 0, "positive_rate": np.nan, "auprc": np.nan, "auroc": np.nan, "sampled_with_replacement": False})
            rows.append(row)
            continue
        if variant == "matched_1to1":
            with_replacement = len(neg) < len(pos)
            used = neg.sample(n=len(pos), random_state=AUPRC_SEED, replace=with_replacement)
        else:
            with_replacement = False
            used = neg
        labels = np.concatenate([np.ones(len(pos)), np.zeros(len(used))])
        scores = 1.0 - np.concatenate([pos[p_col].to_numpy(), used[p_col].to_numpy()])
        precision, recall, _ = precision_recall_curve(labels, scores)
        fpr, tpr, _ = roc_curve(labels, scores)
        row.update({
            "n_ntc_used": int(len(used)), "positive_rate": float(len(pos) / (len(pos) + len(used))),
            "auprc": float(auc(recall, precision)), "auroc": float(auc(fpr, tpr)),
            "sampled_with_replacement": bool(with_replacement),
            "_curves": {"precision": precision, "recall": recall, "fpr": fpr, "tpr": tpr},
        })
        rows.append(row)
    return rows


# --- plots ---------------------------------------------------------------------


def _qq_series(p: np.ndarray):
    tested = np.sort(p[np.isfinite(p)])
    n = tested.size
    if n == 0:
        return None
    expected = (np.arange(1, n + 1) - 0.5) / n
    lower = stats.beta.ppf(0.025, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1)
    upper = stats.beta.ppf(0.975, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1)
    clip = lambda x: -np.log10(np.clip(x, P_FLOOR, 1.0))  # noqa: E731
    return clip(expected), clip(tested), clip(lower), clip(upper)


def qq_plot(panels: list[tuple[str, dict]], path: Path, title: str) -> None:
    """One panel per pair set; each method a series; a 95% band from the first series."""
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4.2), dpi=130, squeeze=False)
    for ax, (name, series) in zip(axes[0], panels):
        band_drawn = False
        top = 1.0
        for label, p in series.items():
            qq = _qq_series(np.asarray(p, dtype=float))
            if qq is None:
                continue
            expected, observed, lower, upper = qq
            if not band_drawn:
                ax.fill_between(expected, lower, upper, color="0.85", lw=0, label="95% band")
                band_drawn = True
            ax.plot(expected, observed, ".", ms=3, label=f"{label} (n={observed.size:,})")
            top = max(top, float(expected.max()), float(observed.max()))
        ax.plot([0, top], [0, top], color="black", lw=0.8, ls="--")
        ax.set_xlabel("expected -log10 p (uniform)")
        ax.set_ylabel(f"observed -log10 p (floored at {P_FLOOR:g})")
        ax.set_title(name)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def histogram_plot(panels: list[tuple[str, dict]], path: Path, title: str) -> None:
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.6), dpi=130, squeeze=False)
    for ax, (name, series) in zip(axes[0], panels):
        for label, p in series.items():
            p = np.asarray(p, dtype=float)
            p = p[np.isfinite(p)]
            if p.size:
                ax.hist(p, bins=20, range=(0, 1), histtype="step", lw=1.2, label=f"{label} (n={p.size:,})")
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("p-value")
        ax.set_ylabel("pairs")
        ax.set_title(name)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def scatter_plot(joined_by_set: list[tuple[str, pd.DataFrame]], path: Path, title: str) -> None:
    if not joined_by_set:
        return
    fig, axes = plt.subplots(1, len(joined_by_set), figsize=(4.2 * len(joined_by_set), 4.2), dpi=130, squeeze=False)
    for ax, (name, joined) in zip(axes[0], joined_by_set):
        both = joined.dropna(subset=["sceptre_p_value", "perturbo_p_value"])
        s = -np.log10(np.clip(both["sceptre_p_value"], P_FLOOR, 1.0))
        t = -np.log10(np.clip(both["perturbo_p_value"], P_FLOOR, 1.0))
        ax.plot(s, t, ".", ms=3, alpha=0.6)
        top = max(1.0, float(s.max()) if len(s) else 1.0, float(t.max()) if len(t) else 1.0)
        ax.plot([0, top], [0, top], color="black", lw=0.8, ls="--")
        ax.set_xlabel("SCEPTRE -log10 p")
        ax.set_ylabel("PerTurbo -log10 p")
        ax.set_title(f"{name} (n={len(both):,})")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def curves_plot(rows: list[dict], path: Path, level: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=130)
    drawn = False
    colors = {"sceptre": "C0", "perturbo": "C1"}
    for row in rows:
        curves = row.get("_curves")
        if not curves:
            continue
        drawn = True
        style = "-" if row.get("negatives") == "matched_1to1" else "--"
        label = f"{row['method']}, {row.get('negatives', '')} negatives"
        color = colors.get(row["method"])
        axes[0].plot(curves["recall"], curves["precision"], lw=1.2, ls=style, color=color, label=f"{label} AUPRC={row['auprc']:.3f}")
        axes[1].plot(curves["fpr"], curves["tpr"], lw=1.2, ls=style, color=color, label=f"{label} AUROC={row['auroc']:.3f}")
    axes[0].set_xlabel("recall"); axes[0].set_ylabel("precision")
    axes[1].set_xlabel("false positive rate"); axes[1].set_ylabel("true positive rate")
    axes[1].plot([0, 1], [0, 1], color="black", lw=0.8, ls="--")
    # The curves hug the top of the precision panel and the top-left of the
    # ROC panel, so the legends sit in the corners they leave empty.
    axes[0].legend(fontsize=8, frameon=False, loc="lower left")
    axes[1].legend(fontsize=8, frameon=False, loc="lower right")
    fig.suptitle(f"Direct targets versus non-targeting pairs, {level} level")
    fig.tight_layout()
    if drawn:
        fig.savefig(path)
    plt.close(fig)


def cells_plot(counts: pd.DataFrame, path: Path, title: str, column: str = "n_cells") -> None:
    roles = [r for r in ("targeting", "ntc") if (counts["role"] == r).any()]
    if not roles:
        return
    fig, ax = plt.subplots(figsize=(5, 4), dpi=130)
    data = [np.log10(counts.loc[counts["role"] == r, column].to_numpy(dtype=float) + 1) for r in roles]
    ax.boxplot(data, tick_labels=[f"{r} (n={len(d):,})" for r, d in zip(roles, data)], showfliers=True)
    ax.set_ylabel(f"log10({column} + 1)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --- cell-count summary ---------------------------------------------------------


def cell_count_summary(per_guide: pd.DataFrame, per_element: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    rows, warnings = [], []
    for level, frame in (("guide", per_guide), ("element", per_element)):
        medians = {}
        for role in ("targeting", "ntc"):
            values = frame.loc[frame["role"] == role, "n_cells"].to_numpy(dtype=float)
            if values.size == 0:
                continue
            q1, q3 = np.percentile(values, [25, 75])
            medians[role] = float(np.median(values))
            rows.append({"level": level, "role": role, "n": int(values.size), "median_cells": medians[role],
                         "iqr_low": float(q1), "iqr_high": float(q3), "warning": ""})
        if "targeting" in medians and "ntc" in medians and medians["targeting"] > 0:
            ratio = medians["ntc"] / medians["targeting"]
            if ratio < 0.5 or ratio > 2.0:
                text = (f"{level} level: the non-targeting median of {medians['ntc']:.0f} cells is {ratio:.2f}x the "
                        f"targeting median of {medians['targeting']:.0f}; cell counts differ beyond 0.5x to 2x.")
                warnings.append(text)
                for row in rows:
                    if row["level"] == level:
                        row["warning"] = text
    return pd.DataFrame(rows), warnings


# --- main ---------------------------------------------------------------------------


def _skip(args) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    design = json.loads(Path(args.summary_in).read_text())
    reason = design.get("reason", "the prepare step skipped the calibration pass")
    (outdir / "ntc_calibration_skipped.txt").write_text(f"NTC calibration skipped: {reason}\n")
    per_guide = pd.read_csv(args.cell_counts_per_guide, sep="\t")
    per_element = pd.read_csv(args.cell_counts_per_element, sep="\t")
    summary, warnings = cell_count_summary(per_guide, per_element)
    summary.to_csv(outdir / "ntc_cell_count_summary.tsv", sep="\t", index=False)
    cells_plot(per_guide, outdir / "ntc_cells_per_guide.png", "cells per guide")
    cells_plot(per_element, outdir / "ntc_cells_per_element.png", "cells per element")
    (outdir / "ntc_calibration_notes.txt").write_text("\n".join([f"NTC calibration skipped: {reason}", *warnings]) + "\n")
    print(f"NTC calibration skipped: {reason}", flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.skipped:
        return _skip(args)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    design = json.loads(Path(args.summary_in).read_text())
    mode = json.loads(Path(args.inference_mode).read_text()) if args.inference_mode else {}
    notes: list[str] = []

    pairs = pd.read_csv(args.pairs, sep="\t", dtype={"guide_id": str, "gene_id": str})
    pairs["ntc_bucket"] = pairs["ntc_bucket"].astype(int)
    pairs["pair_set"] = pairs["pair_set"].astype(str)
    element_pairs = pairs.drop_duplicates(["ntc_bucket", "gene_id"])[["ntc_bucket", "gene_id", "pair_set"]]
    guide_targets = pd.read_csv(args.guide_targets, sep="\t", dtype={"guide_id": str})
    direct = pd.read_csv(args.direct_targets, sep="\t", dtype={"guide_id": str, "gene_id": str})
    ntc_guides = set(guide_targets.loc[guide_targets["role"] == "ntc", "guide_id"].astype(str))
    bucket_of_guide = dict(zip(
        guide_targets.loc[guide_targets["role"] == "ntc", "guide_id"].astype(str),
        guide_targets.loc[guide_targets["role"] == "ntc", "ntc_bucket"],
    ))

    # --- the matched tables, per level -----------------------------------------
    sceptre_pass = {
        "guide": load_method_table(args.sceptre_ntc_per_guide, "sceptre", "guide"),
        "element": load_method_table(args.sceptre_ntc_per_element, "sceptre", "element"),
    }
    perturbo_global = {
        "guide": load_method_table(args.perturbo_global_per_guide, "perturbo", "guide"),
        "element": load_method_table(args.perturbo_global_per_element, "perturbo", "element"),
    }
    sceptre_main = {
        "guide": load_method_table(args.sceptre_local_per_guide, "sceptre", "guide") if args.sceptre_local_per_guide else None,
        "element": load_method_table(args.sceptre_local_per_element, "sceptre", "element") if args.sceptre_local_per_element else None,
    }

    # SCEPTRE's element rows are named ntc_calibration|N; PerTurbo's non-targeting|N.
    for level, frame in sceptre_pass.items():
        key = "guide_id" if level == "guide" else "intended_target_name"
        if level == "element":
            frame["ntc_bucket"] = frame[key].map(_bucket)
        else:
            frame["ntc_bucket"] = frame[key].map(bucket_of_guide)
    for level, frame in perturbo_global.items():
        key = "guide_id" if level == "guide" else "intended_target_name"
        if level == "element":
            frame["is_ntc"] = frame[key].str.startswith(CONTROL_PREFIX)
            frame["ntc_bucket"] = frame[key].where(frame["is_ntc"]).map(_bucket)
        else:
            frame["is_ntc"] = frame[key].isin(ntc_guides)
            frame["ntc_bucket"] = frame[key].map(bucket_of_guide)

    joined: dict[str, pd.DataFrame] = {}
    calibration, agreement, auprc = [], [], []
    for level in LEVELS:
        if level == "guide":
            keys = ["guide_id", "gene_id"]
            left = pairs[["guide_id", "gene_id", "ntc_bucket", "pair_set"]].drop_duplicates(keys)
        else:
            keys = ["ntc_bucket", "gene_id"]
            left = element_pairs
        s = sceptre_pass[level].drop(columns=[c for c in ("intended_target_name",) if c in sceptre_pass[level].columns])
        t = perturbo_global[level]
        t = t[t["is_ntc"]].drop(columns=["is_ntc"] + [c for c in ("intended_target_name",) if c in t.columns])
        table = left.merge(s, on=keys, how="left").merge(t, on=keys, how="left", suffixes=("", "_perturbo"))
        if "ntc_bucket_perturbo" in table.columns:
            table = table.drop(columns=["ntc_bucket_perturbo"])
        table = table.drop(columns=[c for c in table.columns if c.endswith("_perturbo") and c != "ntc_bucket_perturbo"], errors="ignore")
        joined[level] = table

        sets = sorted({name for label in table["pair_set"].unique() for name in str(label).split("+")})
        for pair_set in sets:
            subset = table[table["pair_set"].str.split("+").map(lambda parts, ps=pair_set: ps in parts)]
            for method in METHODS:
                column = f"{method}_p_value"
                if column not in subset.columns:
                    continue
                low = subset["perturbo_crt_low_information"] if (method == "perturbo" and "perturbo_crt_low_information" in subset.columns) else None
                calibration.append(calibration_rows(subset[column].to_numpy(dtype=float), level=level, pair_set=pair_set, method=method, low_information=low))
            agreement.append(agreement_row(subset, level=level, pair_set=pair_set))
        # PerTurbo-only: every NTC row of the global table, all genes.
        all_ntc = perturbo_global[level][perturbo_global[level]["is_ntc"]]
        calibration.append(calibration_rows(all_ntc["perturbo_p_value"].to_numpy(dtype=float), level=level, pair_set="perturbo_global_all", method="perturbo",
                                            low_information=all_ntc["perturbo_crt_low_information"] if "perturbo_crt_low_information" in all_ntc.columns else None))

        # --- AUPRC: direct-target positives against NTC negatives on the same genes
        if level == "guide":
            positives_keys = direct[["guide_id", "gene_id"]]
        else:
            element_of_guide = dict(zip(guide_targets["guide_id"].astype(str), guide_targets["intended_target_name"].astype(str)))
            positives_keys = direct.assign(intended_target_name=direct["guide_id"].map(element_of_guide))[["intended_target_name", "gene_id"]].drop_duplicates()
        for method in METHODS:
            source = sceptre_main[level] if method == "sceptre" else perturbo_global[level]
            if source is None:
                notes.append(f"No main-run SCEPTRE {level} table was supplied; its AUPRC is not reported.")
                continue
            pos = positives_keys.merge(source, on=list(positives_keys.columns), how="inner")
            if method == "sceptre":
                neg = sceptre_pass[level]
            else:
                neg = all_ntc
            auprc.extend(auprc_rows(pos, neg, level=level, method=method, p_col=f"{method}_p_value"))

    # --- write tables ------------------------------------------------------------
    per_guide_counts = pd.read_csv(args.cell_counts_per_guide, sep="\t")
    per_element_counts = pd.read_csv(args.cell_counts_per_element, sep="\t")
    guide_cells = dict(zip(per_guide_counts["guide_id"].astype(str), per_guide_counts["n_cells"]))
    element_cells = {int(b): n for b, n in zip(per_element_counts.loc[per_element_counts["role"] == "ntc", "ntc_bucket"].dropna().astype(int),
                                                   per_element_counts.loc[per_element_counts["role"] == "ntc", "n_cells"].dropna())}
    joined["guide"]["n_cells_guide"] = joined["guide"]["guide_id"].map(guide_cells)
    joined["element"]["n_cells_bucket"] = joined["element"]["ntc_bucket"].map(element_cells)
    joined["guide"].to_csv(outdir / "ntc_calibration_per_guide.tsv.gz", sep="\t", index=False, compression="gzip")
    joined["element"].to_csv(outdir / "ntc_calibration_per_element.tsv.gz", sep="\t", index=False, compression="gzip")

    calibration_frame = pd.DataFrame(calibration)
    agreement_frame = pd.DataFrame(agreement)
    calibration_frame.to_csv(outdir / "ntc_calibration_summary.tsv", sep="\t", index=False)
    agreement_frame.to_csv(outdir / "ntc_method_agreement.tsv", sep="\t", index=False)
    auprc_frame = pd.DataFrame([{k: v for k, v in row.items() if k != "_curves"} for row in auprc])
    auprc_frame.to_csv(outdir / "ntc_aucpr_summary.tsv", sep="\t", index=False)
    cells_summary, cell_warnings = cell_count_summary(per_guide_counts, per_element_counts)
    cells_summary.to_csv(outdir / "ntc_cell_count_summary.tsv", sep="\t", index=False)

    ratios = design.get("self_inclusion_ratio_by_bucket", {}) or {}
    leniency = max((float(v) for v in ratios.values() if v is not None), default=None)
    notes.append(f"Regime: {design.get('regime')} (analysis MOI {design.get('analysis_moi')}); reference pool {design.get('reference_pool_size')} cells.")
    if leniency is not None:
        notes.append(
            f"Self-inclusion: each tested non-targeting bucket is at most {leniency:.1%} of the reference it is "
            "compared against, for both methods; the check is lenient by about that much."
        )
    notes.extend(cell_warnings)
    (outdir / "ntc_calibration_notes.txt").write_text("\n".join(notes) + "\n")
    summary = {
        "status": "ok",
        "design": design,
        "inference_mode": {k: mode.get(k) for k in ("analysis_moi", "sceptre_control_group", "perturbo_crt_pool") if k in mode},
        "calibration": calibration_frame.to_dict(orient="records"),
        "agreement": agreement_frame.to_dict(orient="records"),
        "auprc": auprc_frame.to_dict(orient="records"),
        "cell_counts": cells_summary.to_dict(orient="records"),
        "notes": notes,
    }
    (outdir / "ntc_calibration_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # --- plots -------------------------------------------------------------------------
    for level in LEVELS:
        table = joined[level]
        sets = sorted({name for label in table["pair_set"].unique() for name in str(label).split("+")})
        panels, hist_panels, scatter_panels = [], [], []
        for pair_set in sets:
            subset = table[table["pair_set"].str.split("+").map(lambda parts, ps=pair_set: ps in parts)]
            series = {m: subset[f"{m}_p_value"].to_numpy(dtype=float) for m in METHODS if f"{m}_p_value" in subset.columns}
            panels.append((pair_set, series))
            hist_panels.append((pair_set, series))
            scatter_panels.append((pair_set, subset))
        qq_plot(panels, outdir / f"ntc_calibration_qq_{level}.png", f"Non-targeting calibration, {level} level")
        histogram_plot(hist_panels, outdir / f"ntc_calibration_pvalue_hist_{level}.png", f"Non-targeting p-values, {level} level")
        scatter_plot(scatter_panels, outdir / f"ntc_method_scatter_{level}.png", f"Matched non-targeting pairs, {level} level")
        curves_plot([row for row in auprc if row["level"] == level], outdir / f"ntc_aucpr_curves_{level}.png", level)
    all_ntc_series = {level: perturbo_global[level].loc[perturbo_global[level]["is_ntc"], "perturbo_p_value"].to_numpy(dtype=float) for level in LEVELS}
    qq_plot([(level, {"perturbo, all genes": p}) for level, p in all_ntc_series.items()], outdir / "ntc_perturbo_global_qq.png", "PerTurbo non-targeting rows over all genes")
    cells_plot(per_guide_counts, outdir / "ntc_cells_per_guide.png", "cells per guide")
    cells_plot(per_element_counts, outdir / "ntc_cells_per_element.png", "cells per element")
    print("NTC calibration report written to", outdir, flush=True)
    for note in notes:
        print(note, flush=True)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--skipped", action="store_true", help="the prepare step skipped the pass; write the skip note and the cell-count outputs only")
    parser.add_argument("--sceptre-ntc-per-guide")
    parser.add_argument("--sceptre-ntc-per-element")
    parser.add_argument("--perturbo-global-per-guide")
    parser.add_argument("--perturbo-global-per-element")
    parser.add_argument("--sceptre-local-per-guide", default=None, help="main-run SCEPTRE per-guide table (AUPRC positives)")
    parser.add_argument("--sceptre-local-per-element", default=None)
    parser.add_argument("--pairs")
    parser.add_argument("--guide-targets")
    parser.add_argument("--direct-targets")
    parser.add_argument("--summary-in", required=True)
    parser.add_argument("--inference-mode", default=None)
    parser.add_argument("--cell-counts-per-guide", required=True)
    parser.add_argument("--cell-counts-per-element", required=True)
    parser.add_argument("--outdir", default="calibration")
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
