#!/usr/bin/env python
# %%
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from muon import MuData

sys.path.append("/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb/src")
import utils

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["font.size"] = 14
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["xtick.labelsize"] = 14
mpl.rcParams["ytick.labelsize"] = 14
mpl.rcParams["legend.fontsize"] = 14
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.transparent"] = True

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
path_to_cc_genes = "/oak/stanford/groups/engreitz/Users/opushkar/common_sc"
plate_map_path = os.path.join(core_path, "0.data_preparation/1.input_data/TF_plate_map.xlsx")
ens2symbol_path = "/oak/stanford/groups/engreitz/Users/opushkar/genome/ensembl_to_symbol_dict_v43.npy"
input_path = os.path.join(core_path, "2.qc_gex/0.demux_data/data")
output_path = os.path.join(core_path, "2.qc_gex/1.qc_data")
save_plots = True
reference = "human"

days = ["d0", "d1", "d2", "d3"]
reps = [1, 2]

min_cells_per_gene = 500        # ~0.05% of 1M cells
# min_counts_per_gene = 500        # low floor, mostly redundant
# filter_cells_by_min_genes = False
# min_genes_per_cell = 200

min_counts_per_cell = 500

# mt/ribo filtering is NOT done per library -- it's applied once, after concatenating all
# 8 day/rep tags together, in 2.0.concatenate_adata.py. Per-library QC here is UMI-count and
# gene-count based only (cell calling below + the adaptive upper bounds further down), plus
# doublet removal.
n_mads_common = None
n_mads_dict = None

# Cell calling (BarcodeRanks_Inflection via DropletUtils). This is the only cell-calling step
# in the whole pipeline -- 0.demultiplex_adata.py no longer applies any per-sub filter, it just
# pools raw barcodes across all 13 sub-libraries per tag. Computed at tag scale (~1.4M raw
# barcodes/tag, all 13 subs pooled) since that's where a genuine ambient-vs-cell transition is
# visible; per-sub scale (~110k raw barcodes/sub) does have a real Knee-detectable transition
# too (verified empirically), but Inflection specifically is a near no-op at that scale.
# Using inflection rather than knee: knee is markedly stricter (totals ~524k cells across all
# 8 tags vs the wet lab's ~650k-1M expectation) while inflection totals ~712k, squarely inside
# that range and (after fixing a float-truncation bug in run_barcode_ranks_knee that previously
# destabilized it) stable across all 8 tags.
rscript_kb_bin = "/home/groups/engreitz/Users/tony/anaconda3/envs/kb/bin/Rscript"
barcode_ranks_script = os.path.join(core_path, "src/barcode_ranks_knee.R")
barcode_ranks_tmp_dir = os.path.join(output_path, "barcode_ranks_tmp")

random_seed = 0
# Scrublet's default (0.05) is calibrated for droplet-based capture, where multiplets come
# from Poisson loading of a physical partition. This is split-pool combinatorial (Parse
# Evercode) barcoding, where a "doublet" instead requires two cells to independently collide
# on the same barcode combination across all rounds -- Parse's own species-mixing validation
# reports ~2.3% multiplets, well below the droplet-oriented default. Using 0.025 here as a
# closer starting point; this only shifts Scrublet's score threshold; the underlying
# doublet_score itself is unaffected.
expected_doublet_rate = 0.025
# Per-sub automatic threshold detection is unreliable at this scale (~60-70k cells/sub): several
# subs show no real bimodal separation between observed and simulated-doublet score distributions,
# so the automatic detector can lock onto a spurious minimum near the left edge (~0.02) and call
# ~99% of cells as doublets, versus ~0.3 on subs where it works correctly. Fixing an explicit
# threshold (matching where auto-detection lands on well-behaved subs) avoids that failure mode.
scrublet_threshold = 0.3

labeling_dict = {
    "n_genes_by_counts": "# genes per cell",
    "total_counts": "# UMIs per cell",
    "pct_counts_mt": "% mitochondrial counts",
    "pct_counts_ribo": "% ribosomal counts",
    "S_score": "S score",
    "G2M_score": "G2M score",
    "bci": "Barcode index (bci)",
    "stype": "Round1 barcode type (stype)",
    "doublet_score": "Doublet score",
}

plate_map = pd.read_excel(plate_map_path)
plate_map["sample_key"] = plate_map["Sample"].str.replace("_", "")  # d0rep1, d1rep2, etc.

sub_colors_list = plt.cm.tab20(np.linspace(0, 1, 13))
well_colors_list = plt.cm.tab20(np.linspace(0, 1, 12))

np.random.seed(random_seed)

# %%
os.makedirs(os.path.join(output_path, "qc_filtered_data"), exist_ok=True)
ens2symbol = np.load(ens2symbol_path, allow_pickle=True).item()

# qc_thresholds holds n_mads (not raw values) -- is_outlier() and umi_gene_scatter() both
# already expect n_mads here; adaptive per-tag upper bounds are derived from these below.
qc_thresholds = utils.get_qc_thresholds(n_mads_common=n_mads_common, n_mads_dict=n_mads_dict)
#%%

print("QC n_mads thresholds:")
for k, v in qc_thresholds.items():
    print(f"  {k}: {v}")

# %%
# Each day/rep tag is independent (own h5mu, own Scrublet runs, own output file), so under a
# Slurm array job (submit_1_qc.sh, --array=0-7) each task processes exactly one tag instead of
# looping over all 8 serially. Running the script directly (no SLURM_ARRAY_TASK_ID) falls back
# to the original sequential behavior over all tags.
array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
tag_list = [(d, r) for d in days for r in reps]
if array_task_id is not None:
    day_for_task, rep_for_task = tag_list[int(array_task_id)]
    days_to_run, reps_to_run = [day_for_task], [rep_for_task]
else:
    days_to_run, reps_to_run = days, reps

summary_rows = []

for day in days_to_run:
    for rep in reps_to_run:
        tag = f"{day}_rep{rep}"
        plot_dir = os.path.join(output_path, "plots", "QC", tag)
        os.makedirs(plot_dir, exist_ok=True)

        h5mu_path = os.path.join(input_path, f"{tag}_all_subs_gex_and_guide.h5mu")
        assert os.path.exists(h5mu_path), f"Missing h5mu: {h5mu_path}"
        print(f"\nLoading {tag}")

        mdata = mu.read(h5mu_path)
        gex_adata = mdata["GEX"]
        guide_adata = mdata["guide"]
        gex_adata.X = gex_adata.X.astype(np.float32)

        subs = sorted(gex_adata.obs["sub"].unique(), key=lambda x: int(x))
        sub_colors = {s: sub_colors_list[i] for i, s in enumerate(subs)}

        # well colors: well positions for this day+rep (plate_map Sample matches tag format d0_rep1)
        wells = sorted(plate_map[plate_map["Sample"] == tag]["Well Position"].tolist())
        well_colors = {w: well_colors_list[i] for i, w in enumerate(wells)}

        # Knee plot
        knee_df = pd.DataFrame(
            {"sum": np.array(gex_adata.X.sum(axis=1)).flatten(), "barcodes": gex_adata.obs_names.values}
        )
        knee_df = knee_df.sort_values("sum", ascending=False).reset_index(drop=True)
        knee_df["sum_log"] = np.log1p(knee_df["sum"])

        # Cell calling: BarcodeRanks knee + inflection (DropletUtils), computed at tag scale
        # (see rscript_kb_bin definition above for why this differs from the perturb_pipeline
        # run); inflection is the one actually used for filtering below, see comment there.
        knee_threshold, inflection_threshold = utils.run_barcode_ranks_knee(
            knee_df["sum"].values, barcode_ranks_tmp_dir, tag, rscript_kb_bin, barcode_ranks_script
        )
        print(f"  BarcodeRanks knee={knee_threshold:.0f}, inflection={inflection_threshold:.0f}")

        if save_plots:
            utils.knee_plot(
                knee_df, experiment=tag, save=True, output_path=plot_dir,
                knee_threshold=knee_threshold, inflection_threshold=inflection_threshold,
                used_threshold="inflection",
            )
            knee_df.to_csv(os.path.join(plot_dir, f"{tag}_knee_plot_data.tsv"), sep="\t", index=False)

        # Gene annotation
        gex_adata.var.reset_index(inplace=True)
        gex_adata.var.index = gex_adata.var["gene_id"]
        gex_adata.var["symbol"] = gex_adata.var["gene_id"].str.split(".").str[0].map(ens2symbol)

        if reference == "human":
            mt_prefix, ribo_prefix = "MT-", ("RPS", "RPL")
        else:
            mt_prefix, ribo_prefix = "Mt-", ("Rps", "Rpl")

        gex_adata.var["mt"] = gex_adata.var["symbol"].str.startswith(mt_prefix)
        gex_adata.var["ribo"] = gex_adata.var["symbol"].str.startswith(ribo_prefix)

        sc.pp.calculate_qc_metrics(gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

        # Cell cycle scoring
        if os.path.exists(path_to_cc_genes):
            gex_adata.var_names = gex_adata.var["symbol"].astype(str)
            gex_adata.var_names_make_unique()
            s_genes, g2m_genes = utils.get_cell_cycle_genes(reference, path_to_cc_genes)
            if s_genes is not None and g2m_genes is not None:
                sc.tl.score_genes_cell_cycle(gex_adata, s_genes=s_genes, g2m_genes=g2m_genes)
                if save_plots:
                    utils.cell_cycle_barplot(
                        gex_adata,
                        experiment=tag,
                        output_dir=plot_dir,
                        file_name=f"{tag}_cell_cycle_before_filtering",
                        save=True,
                    )
            gex_adata.var_names = gex_adata.var["gene_id"]

        print(f"Before filtering: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

        # QC plots before filtering
        if save_plots:
            utils.qc_violin_plot(
                gex_adata,
                qc_thresholds=qc_thresholds,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter",
                save=True,
                # pct_counts_mt/pct_counts_ribo dropped: not filtered per library anymore
                # (moved to 2.0.concatenate_adata.py, post-concatenation)
                qc_metrics=["n_genes_by_counts", "total_counts"],
                output_path=plot_dir,
            )
            utils.umi_gene_scatter(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter",
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter_with_thresholds",
                qc_thresholds=qc_thresholds,
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_categorical(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter",
                feature="well",
                category_colors=well_colors,
                qc_thresholds=qc_thresholds,
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_with_feature(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter",
                cbar_label=labeling_dict["bci"],
                feature="bci",
                qc_thresholds=qc_thresholds,
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_with_feature(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_before_filter",
                cbar_label=labeling_dict["stype"],
                feature="stype",
                qc_thresholds=qc_thresholds,
                save=True,
                output_path=plot_dir,
            )

        # Cell calling (BarcodeRanks_Inflection): drop ambient/empty-droplet barcodes below the
        # inflection point before any other filtering. Applied after the before_filter plots
        # above so those plots still show the full raw barcode set for diagnostic purposes.
        n_before_cell_calling = gex_adata.n_obs
        called_barcodes = knee_df.loc[knee_df["sum"] >= inflection_threshold, "barcodes"].values
        gex_adata = gex_adata[called_barcodes].copy()
        n_after_cell_calling = gex_adata.n_obs
        print(f"  cell calling (inflection): {n_before_cell_calling} → {n_after_cell_calling} barcodes")

        # Gene filter moved ahead of Scrublet (below) so doublet detection runs on a much
        # smaller matrix -- Scrublet doesn't need the full ~60k-gene raw feature space.
        n_genes_before_filter = gex_adata.n_vars
        sc.pp.filter_genes(gex_adata, min_cells=min_cells_per_gene)
        sc.pp.filter_genes(gex_adata, min_counts=min_counts_per_gene)
        print(f"  gene filter: {n_genes_before_filter} → {gex_adata.n_vars} genes")

        # Refresh per-cell QC metrics after removing genes -- total_counts/n_genes_by_counts
        # (and their log1p versions, used below for the adaptive caps) would otherwise still
        # reflect the pre-gene-filter matrix.
        sc.pp.calculate_qc_metrics(gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

        # --- Doublet detection (Scrublet), removed before the other QC filters below ---
        # Run per sub-library, not per tag: a full tag is ~90k called cells, which is outside
        # the cell-count range Scrublet's simulation/PCA/kNN pipeline is designed for (typical
        # per-10x-lane usage is 1k-50k cells). This is also the methodologically correct
        # granularity, not just a compute-budget workaround: a split-pool barcode collision
        # only produces a mixed/doublet-like transcriptome if the two colliding physical cells
        # both land in the *same* sub-library aliquot during the post-barcoding physical split
        # -- if they land in different subs, each aliquot sees a pure singlet under that shared
        # barcode. So doublet status is a property of a sub's own barcode index space, and
        # per-sub Scrublet is the more correct granularity, not merely a way to shrink the run.
        # threshold=scrublet_threshold is passed explicitly rather than relying on Scrublet's
        # automatic threshold detection, which is unreliable at this scale (see scrublet_threshold
        # definition above).
        gex_adata.obs["doublet_score"] = np.nan
        gex_adata.obs["predicted_doublet"] = False
        for sub in subs:
            sub_mask = (gex_adata.obs["sub"] == sub).values
            adata_sub = gex_adata[sub_mask].copy()
            sc.pp.scrublet(
                adata_sub,
                expected_doublet_rate=expected_doublet_rate,
                threshold=scrublet_threshold,
                random_state=random_seed,
            )
            gex_adata.obs.loc[sub_mask, "doublet_score"] = adata_sub.obs["doublet_score"].values
            gex_adata.obs.loc[sub_mask, "predicted_doublet"] = adata_sub.obs["predicted_doublet"].values
            if save_plots:
                sc.pl.scrublet_score_distribution(adata_sub, show=False)
                utils.save_fig(plot_dir, f"{tag}_sub{sub}_scrublet_score_distribution")

        n_before_doublets = gex_adata.n_obs
        n_doublets = int(gex_adata.obs["predicted_doublet"].sum())
        if save_plots:
            gex_adata.obs[["sub", "doublet_score", "predicted_doublet"]].to_csv(
                os.path.join(plot_dir, f"{tag}_scrublet_score_distribution_data.tsv"), sep="\t"
            )
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.hist(
                gex_adata.obs.loc[~gex_adata.obs["predicted_doublet"], "doublet_score"],
                bins=50, alpha=0.7, label="singlet (predicted)", color="lightgray",
            )
            ax.hist(
                gex_adata.obs.loc[gex_adata.obs["predicted_doublet"], "doublet_score"],
                bins=50, alpha=0.7, label="doublet (predicted)", color="#D19FC7",
            )
            ax.set_xlabel("Doublet score")
            ax.set_ylabel("# cells")
            ax.legend()
            utils.save_fig(plot_dir, f"{tag}_scrublet_score_distribution_pooled")
        gex_adata = gex_adata[~gex_adata.obs["predicted_doublet"]].copy()
        n_after_doublets = gex_adata.n_obs
        print(f"  doublet filter: {n_before_doublets} → {n_after_doublets} ({n_doublets} predicted doublets removed)")

        # Adaptive gene/UMI-count upper bounds (per-tag, MAD-based via is_outlier -- see
        # qc_thresholds definition above). Only the upper bound is used: the lower bound is
        # already handled by inflection-based cell calling above. A single fixed cap (e.g.
        # 100000 UMIs) doesn't fit all tags equally well since typical depth already varies
        # ~4x across tags (d0 inflection ~9-10k vs d3 ~4-5k).
        n_before = gex_adata.n_obs
        _, _, total_counts_upper = utils.is_outlier(
            gex_adata, "log1p_total_counts", qc_thresholds["log1p_total_counts"], verbose=True
        )
        _, _, n_genes_upper = utils.is_outlier(
            gex_adata, "log1p_n_genes_by_counts", qc_thresholds["log1p_n_genes_by_counts"], verbose=True
        )
        total_counts_upper = np.expm1(total_counts_upper)
        n_genes_upper = np.expm1(n_genes_upper)
        gex_adata = gex_adata[gex_adata.obs["n_genes_by_counts"] <= n_genes_upper].copy()
        gex_adata = gex_adata[gex_adata.obs["total_counts"] <= total_counts_upper].copy()
        n_after_adaptive_caps = gex_adata.n_obs

        if filter_cells_by_min_genes:
            sc.pp.filter_cells(gex_adata, min_genes=min_genes_per_cell)
            sc.pp.filter_cells(gex_adata, min_counts=min_counts_per_cell)

        print(f"  adaptive gene/UMI caps (total_counts<={total_counts_upper:.0f}, "
              f"n_genes<={n_genes_upper:.0f}): {n_before} → {n_after_adaptive_caps}")
        print(f"Final: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

        # QC plots after filtering
        if save_plots:
            sc.pp.calculate_qc_metrics(gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

            utils.umi_gene_scatter(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_after_filter",
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_categorical(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_after_filter",
                feature="well",
                category_colors=well_colors,
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_with_feature(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_after_filter",
                cbar_label=labeling_dict["bci"],
                feature="bci",
                save=True,
                output_path=plot_dir,
            )
            utils.umi_gene_scatter_with_feature(
                gex_adata,
                labeling_dict=labeling_dict,
                experiment=f"{tag}_after_filter",
                cbar_label=labeling_dict["stype"],
                feature="stype",
                save=True,
                output_path=plot_dir,
            )

            if os.path.exists(path_to_cc_genes):
                gex_adata.var_names = gex_adata.var["symbol"].astype(str)
                gex_adata.var_names_make_unique()
                s_genes, g2m_genes = utils.get_cell_cycle_genes(reference, path_to_cc_genes)
                if s_genes is not None and g2m_genes is not None:
                    sc.tl.score_genes_cell_cycle(gex_adata, s_genes=s_genes, g2m_genes=g2m_genes)
                    utils.cell_cycle_barplot(
                        gex_adata,
                        experiment=tag,
                        output_dir=plot_dir,
                        file_name=f"{tag}_cell_cycle_after_filtering",
                        save=True,
                    )
                gex_adata.var_names = gex_adata.var["gene_id"]

        # Sync guide to filtered GEX barcodes
        filtered_barcodes = list(set(gex_adata.obs_names).intersection(guide_adata.obs_names))
        mdata_filtered = MuData({
            "GEX": gex_adata[filtered_barcodes].copy(),
            "guide": guide_adata[filtered_barcodes].copy(),
        })
        out_path = os.path.join(output_path, "qc_filtered_data", f"{tag}_qc_filtered_gex_and_guide.h5mu")
        mdata_filtered.write(out_path)
        print(f"Written: {out_path}")

        row = {
            "tag": tag, "day": day, "rep": rep,
            "knee_threshold": knee_threshold,
            "inflection_threshold": inflection_threshold,
            "n_barcodes_raw": n_before_cell_calling,
            "n_cells_after_cell_calling": n_after_cell_calling,
            "n_cells_before_doublets": n_before_doublets,
            "n_doublets_removed": n_doublets,
            "n_cells_before_adaptive_caps": n_before,
            "total_counts_upper_bound": total_counts_upper,
            "n_genes_upper_bound": n_genes_upper,
            "n_cells_after_adaptive_caps": n_after_adaptive_caps,
            "n_cells_final": gex_adata.n_obs,
            "n_genes_final": gex_adata.n_vars,
        }
        summary_rows.append(row)
        # Written per tag so parallel array tasks never write the combined summary concurrently;
        # 1.1.aggregate_qc_summary.py merges these once all array tasks finish.
        pd.DataFrame([row]).to_csv(
            os.path.join(output_path, f"qc_filtering_summary_{tag}.tsv"), sep="\t", index=False
        )

# %%
if array_task_id is None:
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_path, "qc_filtering_summary.tsv"), sep="\t", index=False)
    print("\nQC filtering summary:")
    print(summary_df.to_string(index=False))
# %%
