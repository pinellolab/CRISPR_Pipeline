#!/usr/bin/env python
# %%
import os
import subprocess
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

sys.path.append("/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb/src")
import utils

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
qc_filtered_path = os.path.join(core_path, "2.qc_gex/1.qc_data/qc_filtered_data")
guide_reference_path = os.path.join(core_path, "0.data_preparation/1.input_data/TF_tri_guides_qc_reference.txt")
ens2symbol_path = "/oak/stanford/groups/engreitz/Users/opushkar/genome/ensembl_to_symbol_dict_v43.npy"
# Fallback symbol->Ensembl-ID source for targets not in ens2symbol_path (mostly old/renamed
# HGNC symbols, e.g. AARS->AARS1, ARNTL->BMAL1/renamed) -- has its own "gene" (old symbol) ->
# "ensembl_id" mapping that resolves ~54/79 of the otherwise-unmapped targets in this guide
# library (verified 2026-09-04). The remaining ~25 are guide-strength suffixes
# (CD151_strong/_weak), individual-guide-naming artifacts (POLR1D_sgRNA_B), or lncRNA/pseudogene
# IDs not present in either reference -- not real gene symbols, not fixable this way.
weissman_guides_path = os.path.join(core_path, "0.data_preparation/1.input_data/weissman_guides_with_coordinates.tsv")
output_path = os.path.join(core_path, "2.qc_gex/1.qc_data/guide_assignment_data")
tmp_root = os.path.join(output_path, "tmp")

days = ["d0", "d1", "d2", "d3"]
reps = [1, 2]

# Run as its own pass, independent of transcriptome QC/concatenation (see project memory):
# guide assignment is computationally its own thing (per-guide mixture-model fitting), and
# keeping it separate means 2.0.concatenate_adata.py only ever loads the small guide
# AnnData produced here, not a second copy of the large GEX matrix.
#
# Pooled per differentiation DAY (both reps together), not per day-rep tag (2026-09-04
# decision): sceptre's mixture-model EM fit for each guide uses all cells' UMI counts for that
# guide, so pooling ~2x more cells per guide (both reps) gives the EM fit more power to
# separate background from true perturbation, especially for sparsely-represented guides.
rscript_bin = "/home/groups/engreitz/Users/opushkar/.conda/envs/op_pah/bin/Rscript"
assign_grnas_script = os.path.join(core_path, "src/assign_grnas_sceptre.R")

# Multiplicity_of_infection = 'high' in CRISPR_pipeline/conf/tf_2400_tri_sceptre.config for
# this dataset; probability_threshold/n_em_rep left at sceptre's own defaults there too.
moi = "high"
probability_threshold = "default"
n_em_rep = "default"

# Extra covariates for sceptre's mixture model / discovery analysis, on top of sceptre's own
# automatically-derived covariates (response/grna n_umis, n_nonzero, etc.) -- user-requested
# (2026-09-04): sub-library, replicate, %ribosomal counts, %mt counts.
covariate_cols = ["sub", "rep", "pct_counts_ribo", "pct_counts_mt"]

# Discovery analysis: on-target only (each real TF target vs its own gene, ~2400 pairs for this
# dataset -- a knockdown-efficacy check), not genome-wide trans discovery (2026-09-04 decision:
# that would be tens of millions of tests, out of scope for now).
run_discovery = True

# Per-gRNA assignment diagnostic, gRNA count distribution, calibration check, and (if
# run_discovery) discovery analysis plots/tables -- write_outputs_to_directory() in
# assign_grnas_sceptre.R saves the full standard sceptre bundle, one set per day.
save_plots = True
n_grnas_to_plot = 9
n_calibration_pairs = 5000
calibration_group_size = 3  # matches this dataset's "TF tri-guide" (3 guides/target) design
plot_path = os.path.join(core_path, "2.qc_gex/1.qc_data/guide_assignment_data/plots")

os.makedirs(output_path, exist_ok=True)
os.makedirs(tmp_root, exist_ok=True)
if save_plots:
    os.makedirs(plot_path, exist_ok=True)


# %%
# build_intended_target_name/build_discovery_pairs moved to src/utils.py so 2.0.concatenate_adata.py
# and 2.1.integrate_scvi.py can reuse the same target-gene resolution (e.g. to force perturbation
# target genes into the HVG set) without duplicating this logic.
def write_sceptre_inputs(gex_adata, guide_adata, intended_target_name, discovery_pairs, tmp_dir):
    os.makedirs(tmp_dir, exist_ok=True)
    cell_barcodes = gex_adata.obs_names.values
    assert (guide_adata.obs_names.values == cell_barcodes).all(), "GEX/guide barcode order mismatch"

    response_matrix = sp.csr_matrix(gex_adata.X).T  # genes x cells
    sio.mmwrite(os.path.join(tmp_dir, "response_matrix.mtx"), response_matrix)
    with open(os.path.join(tmp_dir, "gene_ids.txt"), "w") as f:
        f.write("\n".join(gex_adata.var_names))

    grna_matrix = sp.csr_matrix(guide_adata.X).T  # guides x cells
    sio.mmwrite(os.path.join(tmp_dir, "grna_matrix.mtx"), grna_matrix)
    with open(os.path.join(tmp_dir, "grna_ids.txt"), "w") as f:
        f.write("\n".join(guide_adata.var_names))

    pd.DataFrame({"grna_id": guide_adata.var_names, "grna_target": intended_target_name}).to_csv(
        os.path.join(tmp_dir, "grna_target.tsv"), sep="\t", index=False
    )
    with open(os.path.join(tmp_dir, "cell_barcodes.txt"), "w") as f:
        f.write("\n".join(cell_barcodes))

    gex_adata.obs[covariate_cols].to_csv(
        os.path.join(tmp_dir, "extra_covariates.tsv"), sep="\t", index=False
    )
    discovery_pairs.to_csv(os.path.join(tmp_dir, "discovery_pairs.tsv"), sep="\t", index=False)


def run_assign_grnas_sceptre(
    tmp_dir, moi, probability_threshold, n_em_rep, save_plots=False, plot_dir=None,
    n_grnas_to_plot=9, n_calibration_pairs=5000, calibration_group_size=3, run_discovery=False,
):
    output_mtx = os.path.join(tmp_dir, "guide_assignment.mtx")
    cmd = [
        rscript_bin, assign_grnas_script,
        "--response_matrix", os.path.join(tmp_dir, "response_matrix.mtx"),
        "--gene_ids", os.path.join(tmp_dir, "gene_ids.txt"),
        "--grna_matrix", os.path.join(tmp_dir, "grna_matrix.mtx"),
        "--grna_ids", os.path.join(tmp_dir, "grna_ids.txt"),
        "--grna_target", os.path.join(tmp_dir, "grna_target.tsv"),
        "--cell_barcodes", os.path.join(tmp_dir, "cell_barcodes.txt"),
        "--extra_covariates", os.path.join(tmp_dir, "extra_covariates.tsv"),
        "--discovery_pairs", os.path.join(tmp_dir, "discovery_pairs.tsv"),
        "--moi", moi, "--output_mtx", output_mtx,
        "--probability_threshold", probability_threshold, "--n_em_rep", n_em_rep,
    ]
    if save_plots:
        cmd += [
            "--save_plots", str(save_plots).upper(), "--plot_dir", plot_dir,
            "--n_grnas_to_plot", str(n_grnas_to_plot),
            "--n_calibration_pairs", str(n_calibration_pairs),
            "--calibration_group_size", str(calibration_group_size),
            "--run_discovery", str(run_discovery).upper(),
        ]
    subprocess.run(cmd, check=True)
    assignment = sio.mmread(output_mtx).T.tocsr()  # back to cells x guides
    return assignment


# %%
# One Slurm array task = one differentiation day (both reps pooled); falls back to the full
# sequential loop when run directly (no SLURM_ARRAY_TASK_ID).
array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
days_to_run = [days[int(array_task_id)]] if array_task_id is not None else days

summary_rows = []

for day in days_to_run:
    tmp_dir = os.path.join(tmp_root, day)
    print(f"\nGuide assignment: {day} (reps {reps} pooled)")

    gex_list, guide_list = [], []
    for rep in reps:
        tag = f"{day}_rep{rep}"
        h5mu_path = os.path.join(qc_filtered_path, f"{tag}_qc_filtered_gex_and_guide.h5mu")
        assert os.path.exists(h5mu_path), f"Missing h5mu: {h5mu_path}"
        mdata = mu.read(h5mu_path)
        gex_list.append(mdata["GEX"])
        guide_list.append(mdata["guide"])

    gex_adata = ad.concat(gex_list, join="outer", fill_value=0)
    guide_adata = ad.concat(guide_list, join="outer", fill_value=0, merge="same")
    del gex_list, guide_list
    print(f"  pooled: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes "
          f"(rep composition: {gex_adata.obs['rep'].value_counts().to_dict()})")

    intended_target_name = utils.build_intended_target_name(guide_adata.var_names, guide_reference_path)
    discovery_pairs = utils.build_discovery_pairs(
        intended_target_name, gex_adata.var_names, ens2symbol_path, weissman_guides_path
    )
    print(f"  discovery pairs (on-target): {len(discovery_pairs)}")

    write_sceptre_inputs(gex_adata, guide_adata, intended_target_name, discovery_pairs, tmp_dir)
    tag_plot_dir = os.path.join(plot_path, day) if save_plots else None
    guide_assignment = run_assign_grnas_sceptre(
        tmp_dir, moi, probability_threshold, n_em_rep,
        save_plots=save_plots, plot_dir=tag_plot_dir, n_grnas_to_plot=n_grnas_to_plot,
        n_calibration_pairs=n_calibration_pairs, calibration_group_size=calibration_group_size,
        run_discovery=run_discovery,
    )

    guide_adata.layers["guide_assignment"] = guide_assignment
    guide_adata.var["intended_target_name"] = intended_target_name

    n_cells_assigned = int((np.asarray(guide_assignment.sum(axis=1)).flatten() > 0).sum())
    print(f"  cells with >=1 assigned guide: {n_cells_assigned} / {guide_adata.n_obs}")

    # Lightweight output: guide AnnData only (small, cells x ~17k guides) -- 2.0 merges this
    # back onto each rep's GEX data (subsetting by barcode) rather than duplicating the large
    # GEX matrix here.
    out_path = os.path.join(output_path, f"{day}_guide_assignment.h5ad")
    guide_adata.write(out_path)
    print(f"Written: {out_path}")

    summary_rows.append({
        "day": day,
        "n_cells": guide_adata.n_obs,
        "n_guides": guide_adata.n_vars,
        "n_cells_assigned": n_cells_assigned,
        "n_discovery_pairs": len(discovery_pairs),
    })
    pd.DataFrame([summary_rows[-1]]).to_csv(
        os.path.join(output_path, f"guide_assignment_summary_{day}.tsv"), sep="\t", index=False
    )

# %%
if array_task_id is None:
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_path, "guide_assignment_summary.tsv"), sep="\t", index=False)
    print("\nGuide assignment summary:")
    print(summary_df.to_string(index=False))
# %%
