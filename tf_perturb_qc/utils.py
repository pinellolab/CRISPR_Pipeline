import os
import subprocess

import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from matplotlib.collections import PathCollection

matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["font.size"] = 14
matplotlib.rcParams["axes.labelsize"] = 14
matplotlib.rcParams["axes.titlesize"] = 14
matplotlib.rcParams["xtick.labelsize"] = 14
matplotlib.rcParams["ytick.labelsize"] = 14
matplotlib.rcParams["legend.fontsize"] = 14
matplotlib.rcParams["figure.dpi"] = 100
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["savefig.bbox"] = "tight"
matplotlib.rcParams["savefig.transparent"] = True


def is_outlier(adata, metric, nmads, verbose=False):
    M = adata.obs[metric]
    median = np.median(M)
    mad = scipy.stats.median_abs_deviation(M)
    lower_bound = median - nmads * mad
    upper_bound = median + nmads * mad

    if "mt" in metric or "ribo" in metric:
        lower_bound = 0

    outlier = (M < lower_bound) | (M > upper_bound)

    if verbose:
        print(f"{metric}: lower={lower_bound:.2f}, upper={upper_bound:.2f}")
        return outlier, lower_bound, upper_bound
    return outlier


def get_qc_thresholds(n_mads_common=None, n_mads_dict=None):
    default_n_mads = {
        "log1p_total_counts": 5,
        "log1p_n_genes_by_counts": 5,
        # "pct_counts_mt": 5,
        # "pct_counts_ribo": 5,
    }
    if n_mads_common is not None:
        thresholds = {k: n_mads_common for k in default_n_mads}
        print(f"Using common n_mads = {n_mads_common}")
    else:
        thresholds = default_n_mads.copy()
        print("Using default n_mads thresholds")
    if n_mads_dict is not None:
        thresholds.update(n_mads_dict)
        print(f"Metric-specific overrides: {n_mads_dict}")
    return thresholds


def get_cell_cycle_genes(reference, path_to_cc_genes):
    prefix_map = {"human": "hs", "mouse": "mm"}
    prefix = prefix_map.get(reference)
    if prefix is None:
        print(f"No cell cycle genes for reference '{reference}'. Skipping.")
        return None, None
    s_genes = pd.read_csv(
        os.path.join(path_to_cc_genes, f"{prefix}_cell_cycle_s_genes.txt"),
        header=None,
        names=["s_genes"],
    ).values.flatten()
    g2m_genes = pd.read_csv(
        os.path.join(path_to_cc_genes, f"{prefix}_cell_cycle_g2m_genes.txt"),
        header=None,
        names=["g2m_genes"],
    ).values.flatten()
    return s_genes, g2m_genes


def save_fig(output_dir, filename):
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(output_dir, f"{filename}.{ext}"))
    plt.close()


def cell_cycle_barplot(
    adata,
    file_name,
    experiment,
    output_dir=None,
    save=False,
    cc_cats=["G1", "S", "G2M"],
    colors=["#15616D", "#FFECD1", "#FF7D00"],
):
    counts = adata.obs.phase.value_counts()
    data = np.array([[counts.get(cat, 0) / adata.obs.shape[0] * 100 for cat in cc_cats]])

    fig, ax = plt.subplots(1, 1, figsize=(1, 5))
    bottom = 0
    for i, cat in enumerate(cc_cats):
        ax.bar(0, data[0, i], bottom=bottom, color=colors[i], label=cat, width=0.15)
        bottom += data[0, i]

    ax.text(-0.07, 102.5, f"n={adata.obs.shape[0]}")
    ax.set_ylabel("% cells")
    ax.set_ylim(0, 100)
    ax.set_xticks([0])
    ax.set_xticklabels([experiment])
    ax.legend(bbox_to_anchor=(1.05, 0.5))
    if save:
        os.makedirs(output_dir, exist_ok=True)
        save_fig(output_dir, file_name)


def knee_plot(
    knee_df, experiment, save=False, output_path=None,
    knee_threshold=None, inflection_threshold=None, used_threshold=None,
):
    """
    used_threshold: which of "knee"/"inflection" is actually used for cell-calling filtering
    downstream (drawn solid instead of dashed, and called out in the legend) -- both are shown
    since 1.qc_adata.py computes both but only filters on one, and it's easy for the plotted
    reference line to silently drift out of sync with whichever one filtering actually uses.
    """
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.plot(
        np.log1p(knee_df.index),
        knee_df["sum_log"],
        marker="o",
        linestyle="-",
        markersize=2,
        color="gray",
        alpha=0.5,
        rasterized=True,
    )
    ax.set_xlabel("log1p( Barcode rank )")
    ax.set_ylabel("log1p( UMI counts )")
    ax.set_title(experiment)
    thresholds = {"knee": (knee_threshold, "#D19FC7"), "inflection": (inflection_threshold, "#6699CC")}
    for name, (value, color) in thresholds.items():
        if value is None:
            continue
        is_used = name == used_threshold
        ax.axhline(np.log1p(value), color=color, linestyle="-" if is_used else "--", linewidth=1.5 if is_used else 1)
        label = f" {name}={int(value)}" + (" (used)" if is_used else "")
        ax.text(ax.get_xlim()[1], np.log1p(value), label, color=color, va="bottom", ha="right")
    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"{experiment}_knee_plot")


def run_barcode_ranks_knee(total_counts, tmp_dir, tag, rscript_bin, r_script_path):
    """
    DropletUtils::barcodeRanks() knee/inflection cell calling, run via Rscript on a
    per-barcode total-UMI-count vector (see src/barcode_ranks_knee.R). Must be called
    at tag scale (all sub-libraries pooled, ~800k-900k barcodes) -- at per-sub scale
    (~65-70k barcodes) no real knee exists (see tf_perturb project memory).
    """
    os.makedirs(tmp_dir, exist_ok=True)
    counts_file = os.path.join(tmp_dir, f"{tag}_total_counts.txt")
    output_file = os.path.join(tmp_dir, f"{tag}_barcode_ranks.tsv")
    # fmt="%d" would truncate (not round) fractional totals, e.g. 157.99998 -> 157, creating
    # artificial ties in the low-count tail that destabilize barcodeRanks' inflection detection
    np.savetxt(counts_file, np.asarray(total_counts), fmt="%.6f")
    subprocess.run(
        [rscript_bin, r_script_path, "--counts_file", counts_file, "--output_file", output_file],
        check=True,
    )
    result = pd.read_csv(output_file, sep="\t")
    return float(result["knee"].iloc[0]), float(result["inflection"].iloc[0])


def combined_knee_plot(
    pool_data_dict,
    label,
    save=False,
    count_key="total_counts",
    output_path=None,
):
    """
    Combined knee plot across subpools, colored by sub.

    pool_data_dict: {sub_key: (adata_all_common, valid_cell_barcodes_list)}
    """
    subs = list(pool_data_dict.keys())
    colors = plt.cm.tab20(np.linspace(0, 1, len(subs)))
    sub_colors = {s: c for s, c in zip(subs, colors)}

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    total_cells = 0
    total_barcodes = 0

    for sub, (adata, valid_barcodes) in pool_data_dict.items():
        if count_key not in adata.obs.columns:
            adata.obs[count_key] = np.array(adata.X.sum(axis=1)).flatten()

        sorted_counts = np.sort(adata.obs[count_key].values)[::-1]
        valid_set = set(valid_barcodes)
        is_cell = np.array([bc in valid_set for bc in adata.obs_names])
        sorted_idx = np.argsort(adata.obs[count_key].values)[::-1]
        is_cell_sorted = is_cell[sorted_idx]
        barcode_ranks = np.arange(1, len(sorted_counts) + 1)

        ax.scatter(
            barcode_ranks[~is_cell_sorted],
            np.log1p(sorted_counts[~is_cell_sorted]),
            c="lightgray",
            s=0.3,
            alpha=0.3,
            rasterized=True,
        )
        ax.scatter(
            barcode_ranks[is_cell_sorted],
            np.log1p(sorted_counts[is_cell_sorted]),
            c=[sub_colors[sub]],
            s=0.8,
            alpha=0.5,
            rasterized=True,
        )
        total_cells += is_cell.sum()
        total_barcodes += len(adata)

    legend_elements = [matplotlib.patches.Patch(facecolor="lightgray", label="Background", alpha=0.5)]
    for sub in subs:
        legend_elements.append(matplotlib.patches.Patch(facecolor=sub_colors[sub], label=f"sub{sub}"))

    ax.set_xscale("log")
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("log1p( UMI counts )")
    ax.set_title(f"Knee plot, {label}")
    ax.legend(handles=legend_elements, frameon=True, loc="upper right", fontsize=7)
    ax.text(
        0.02, 0.02,
        f"#barcodes: {total_barcodes}\n#valid cells: {total_cells}",
        transform=ax.transAxes,
        verticalalignment="bottom",
    )

    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"knee_plot_{label}_all_subs")


def qc_violin_plot(
    adata,
    qc_thresholds,
    labeling_dict,
    experiment,
    save=False,
    show_qc_thresholds=True,
    qc_metrics=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    output_path=None,
):
    qc_metric_map = {
        "n_genes_by_counts": "log1p_n_genes_by_counts",
        "total_counts": "log1p_total_counts",
        "pct_counts_mt": "pct_counts_mt",
        "pct_counts_ribo": "pct_counts_ribo",
    }

    fig, axes = plt.subplots(1, len(qc_metrics), figsize=(len(qc_metrics) * 2, 5))
    axes = np.atleast_1d(axes)

    for i, (ax, qc_metric) in enumerate(zip(axes, qc_metrics)):
        threshold_metric = qc_metric_map[qc_metric]
        n_mads = qc_thresholds[threshold_metric]

        # is_outlier forces lower_bound=0 for mt/ribo metrics internally; for total_counts/
        # n_genes_by_counts the bounds come back in log1p space and need expm1 to compare
        # against the raw-scale violin plotted below.
        _, low, high = is_outlier(adata, threshold_metric, n_mads, verbose=True)
        if qc_metric in ["total_counts", "n_genes_by_counts"]:
            low, high = np.expm1(low), np.expm1(high)
        threshold_str = f"n_mads={n_mads}"
        low = max(0, low)

        sc.pl.violin(adata, qc_metric, jitter=0.2, ax=ax, color="gray", show=False, size=2, alpha=0.4)
        for coll in ax.collections:
            if isinstance(coll, PathCollection):
                coll.set_facecolor("gray")
                coll.set_edgecolor("none")
                coll.set_alpha(0.4)
                coll.set_zorder(2)
                coll.set_rasterized(True)

        data = adata.obs[qc_metric].dropna()
        median_val = np.median(data)
        ax.hlines(median_val, -0.075, 0.075, color="black", lw=1.5, zorder=4)

        str_ = f"{experiment}, n={adata.n_obs}" if i == 0 else ""
        ax.set_title(
            f"{str_}\n{threshold_str}\nthreshold={high:.2f}\nmedian={median_val:.2f}",
            fontsize=9,
        )
        if show_qc_thresholds:
            ax.axhline(low, ls="--", lw=0.5, color="red", zorder=3)
            ax.axhline(high, ls="--", lw=0.5, color="red", zorder=3)
        ax.set_ylabel(labeling_dict.get(qc_metric, qc_metric))
        ax.set_xticklabels([])
        ax.set_xlabel("")

    plt.tight_layout()
    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"{experiment}_qc_violin")


def umi_gene_scatter(
    adata,
    labeling_dict,
    experiment,
    save=False,
    x_label="total_counts",
    y_label="n_genes_by_counts",
    color_by="pct_counts_mt",
    qc_thresholds=None,
    output_path=None,
):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    x, y, c = adata.obs[x_label], adata.obs[y_label], adata.obs[color_by]
    mask = ~(x.isna() | y.isna() | c.isna())
    x, y, c = x[mask], y[mask], c[mask]

    scatter = ax.scatter(x, y, c=c, cmap="Blues", s=5, alpha=0.4, linewidth=0, rasterized=True)

    if qc_thresholds:
        qc_metric_map = {
            "n_genes_by_counts": "log1p_n_genes_by_counts",
            "total_counts": "log1p_total_counts",
        }
        for qc_metric in ["total_counts", "n_genes_by_counts"]:
            _, low, high = is_outlier(adata, qc_metric_map[qc_metric], qc_thresholds[qc_metric_map[qc_metric]], verbose=True)
            low, high = np.expm1(low), np.expm1(high)
            low = max(0, low)
            if qc_metric == "total_counts":
                ax.axvline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axvline(high, ls="--", lw=0.5, color="gray", zorder=3)
            else:
                ax.axhline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axhline(high, ls="--", lw=0.5, color="gray", zorder=3)

    ax.set_xlabel(labeling_dict.get(x_label, x_label))
    ax.set_ylabel(labeling_dict.get(y_label, y_label))
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(labeling_dict.get(color_by, color_by))
    ax.set_title(f"{experiment}, n={adata.n_obs}")
    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"{experiment}_umi_vs_genes_scatter")


def umi_gene_scatter_categorical(
    adata,
    labeling_dict,
    experiment,
    feature,
    category_colors,
    save=False,
    x_label="total_counts",
    y_label="n_genes_by_counts",
    qc_thresholds=None,
    output_path=None,
):
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    x, y, c = adata.obs[x_label], adata.obs[y_label], adata.obs[feature]
    mask = ~(x.isna() | y.isna() | c.isna())
    x, y, c = x[mask], y[mask], c[mask]

    for cat, color in category_colors.items():
        cat_mask = c == cat
        ax.scatter(
            x[cat_mask], y[cat_mask],
            c=color, s=5, alpha=0.4, linewidth=0, rasterized=True,
            label=f"{cat} (n={int(cat_mask.sum())})",
        )

    if qc_thresholds:
        qc_metric_map = {
            "n_genes_by_counts": "log1p_n_genes_by_counts",
            "total_counts": "log1p_total_counts",
        }
        for qc_metric in ["total_counts", "n_genes_by_counts"]:
            _, low, high = is_outlier(adata, qc_metric_map[qc_metric], qc_thresholds[qc_metric_map[qc_metric]], verbose=True)
            low, high = np.expm1(low), np.expm1(high)
            low = max(0, low)
            if qc_metric == "total_counts":
                ax.axvline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axvline(high, ls="--", lw=0.5, color="gray", zorder=3)
            else:
                ax.axhline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axhline(high, ls="--", lw=0.5, color="gray", zorder=3)

    ax.set_xlabel(labeling_dict.get(x_label, x_label))
    ax.set_ylabel(labeling_dict.get(y_label, y_label))
    ax.legend(frameon=True, fontsize=7, markerscale=2)
    ax.set_title(f"{experiment}, n={adata.n_obs}")
    plt.tight_layout()
    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"{experiment}_{feature}_color_scatter")


def replicate_composition_test(adata, cluster_key, batch_key="rep"):
    """
    Chi-square test of independence between cluster assignment and a batch/replicate label,
    with per-cluster standardized residuals (chi2 contribution) and BH-FDR correction across
    clusters. Use to flag clusters whose replicate composition deviates from the dataset-wide
    proportion, i.e. clusters that may be driven by a batch effect rather than biology.
    Returns (per_cluster_df, overall_chi2, overall_p, dof).
    """
    from scipy.stats import chi2 as chi2_dist
    from scipy.stats import chi2_contingency
    from statsmodels.stats.multitest import multipletests

    contingency = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key])
    overall_chi2, overall_p, dof, expected = chi2_contingency(contingency)

    residuals = (contingency.values - expected) / np.sqrt(expected)
    per_cluster_chi2 = (residuals**2).sum(axis=1)
    per_cluster_dof = contingency.shape[1] - 1
    per_cluster_p = chi2_dist.sf(per_cluster_chi2, per_cluster_dof)
    _, per_cluster_fdr, _, _ = multipletests(per_cluster_p, method="fdr_bh")

    result = contingency.copy()
    result.columns = [f"n_{c}" for c in result.columns]
    result["chi2"] = per_cluster_chi2
    result["p_raw"] = per_cluster_p
    result["p_fdr_bh"] = per_cluster_fdr
    result.index.name = cluster_key
    return result.reset_index(), overall_chi2, overall_p, dof


def umi_gene_scatter_with_feature(
    adata,
    labeling_dict,
    experiment,
    cbar_label,
    feature,
    save=False,
    x_label="total_counts",
    y_label="n_genes_by_counts",
    qc_thresholds=None,
    output_path=None,
):
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    x, y, raw_feature = adata.obs[x_label], adata.obs[y_label], adata.obs[feature]
    mask = ~(x.isna() | y.isna() | raw_feature.isna())
    x, y, raw_feature = x[mask], y[mask], raw_feature[mask]

    if pd.api.types.is_numeric_dtype(raw_feature):
        c = raw_feature
        categories = None
    else:
        categories = sorted(raw_feature.unique())
        code_map = {cat: i for i, cat in enumerate(categories)}
        c = raw_feature.map(code_map)

    scatter = ax.scatter(x, y, c=c, cmap="tab20", s=10, alpha=0.5, linewidth=0, rasterized=True)

    if qc_thresholds:
        qc_metric_map = {
            "n_genes_by_counts": "log1p_n_genes_by_counts",
            "total_counts": "log1p_total_counts",
        }
        for qc_metric in ["total_counts", "n_genes_by_counts"]:
            _, low, high = is_outlier(adata, qc_metric_map[qc_metric], qc_thresholds[qc_metric_map[qc_metric]], verbose=True)
            low, high = np.expm1(low), np.expm1(high)
            low = max(0, low)
            if qc_metric == "total_counts":
                ax.axvline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axvline(high, ls="--", lw=0.5, color="gray", zorder=3)
            else:
                ax.axhline(low, ls="--", lw=0.5, color="gray", zorder=3)
                ax.axhline(high, ls="--", lw=0.5, color="gray", zorder=3)

    ax.set_xlabel(labeling_dict.get(x_label, x_label))
    ax.set_ylabel(labeling_dict.get(y_label, y_label))
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(cbar_label)
    if categories is not None:
        cbar.set_ticks(range(len(categories)))
        cbar.set_ticklabels(categories)
    ax.set_title(f"{experiment}, n={adata.n_obs}")
    plt.tight_layout()
    if save and output_path:
        os.makedirs(output_path, exist_ok=True)
        save_fig(output_path, f"{experiment}_{feature}_colorbar_scatter")


def build_intended_target_name(guide_ids, guide_reference_path):
    """Map each guide ID to its intended target gene symbol (or "non-targeting")."""
    ref = pd.read_csv(guide_reference_path, sep="\t")
    missing = set(guide_ids) - set(ref["ID"])
    assert not missing, f"{len(missing)} guide IDs not found in reference: {list(missing)[:5]}"
    # Strip stray whitespace (e.g. "AARS " vs "AARS" both present in the reference file) --
    # otherwise these would form two separate sceptre grna_target groups for what should be one
    # target, splitting its guide set and weakening the mixture-model fit and discovery test.
    target_map = dict(zip(ref["ID"], ref["gene"].str.strip()))
    targets = pd.Series(guide_ids).map(target_map)
    # sceptre's reserved keyword for controls is the literal string "non-targeting"; the
    # reference file instead gives each non-targeting guide its own pseudo-gene name
    # (e.g. "non-targeting_00642"), so collapse all of those to the literal keyword.
    targets = targets.where(~targets.str.startswith("non-targeting"), "non-targeting")
    return targets.values


def resolve_target_gene_ids(targets, gene_ids, ens2symbol_path, weissman_guides_path):
    """
    Map target gene symbols to their Ensembl gene ID present in gene_ids (versioned), for
    building sceptre discovery pairs or for forcing perturbation-target genes into an HVG set
    regardless of variance rank. ens2symbol_path (gencode v43) is the primary source; the
    weissman guide library's own gene->ensembl_id mapping is a fallback for old/renamed HGNC
    symbols not in ens2symbol_path (recovers ~54/79 otherwise-unmapped targets in this guide
    library, verified 2026-09-04) -- neither source is a strict superset of the other (each
    covers some pseudogenes/lncRNAs/renamed symbols the other doesn't), so both are checked.
    Returns {target: resolved_gene_id} for successfully resolved targets only.
    """
    ens2symbol = np.load(ens2symbol_path, allow_pickle=True).item()
    unversioned_to_versioned = {g.split(".")[0]: g for g in gene_ids}
    symbol_to_versioned = {
        ens2symbol[unv]: versioned
        for unv, versioned in unversioned_to_versioned.items()
        if unv in ens2symbol
    }
    weissman = pd.read_csv(weissman_guides_path, sep="\t", usecols=["gene", "ensembl_id"], low_memory=False)
    weissman_symbol_to_unversioned = (
        weissman.dropna(subset=["gene", "ensembl_id"]).drop_duplicates("gene").set_index("gene")["ensembl_id"]
    )

    def resolve(target):
        target = target.strip()
        if target in symbol_to_versioned:
            return symbol_to_versioned[target]
        unv = weissman_symbol_to_unversioned.get(target)
        if unv is not None and unv in unversioned_to_versioned:
            return unversioned_to_versioned[unv]
        return None

    resolved = {}
    for t in sorted(set(targets)):
        rid = resolve(t)
        if rid is not None:
            resolved[t] = rid
    return resolved


def build_discovery_pairs(intended_target_name, gene_ids, ens2symbol_path, weissman_guides_path):
    """On-target pairs: each real TF target vs its own gene's Ensembl ID in gene_ids."""
    targets = sorted(set(intended_target_name) - {"non-targeting"})
    resolved = resolve_target_gene_ids(targets, gene_ids, ens2symbol_path, weissman_guides_path)
    n_unmapped = len(targets) - len(resolved)
    if n_unmapped:
        print(f"  {n_unmapped}/{len(targets)} targets have no matching gene in the expression matrix, skipped")
    rows = [{"grna_target": t, "response_id": rid} for t, rid in resolved.items()]
    return pd.DataFrame(rows)
