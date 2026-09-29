#!/usr/bin/env python
# %%
"""
Step 2.2: simple concatenation + joint normalization, as a standalone "no batch correction"
alternative to 2.1.integrate_scvi.py -- reads the same mt-filtered GEX+guide checkpoints from
2.0.prepare_adata.py, so no load/concat/guide-filter work is repeated a third time.

Normalization: plain sc.pp.normalize_total(target_sum=None, i.e. median depth) + sc.pp.log1p(),
applied once, globally across all cells for a given mt threshold -- deliberately the simplest
standard scanpy recipe (no per-batch/covariate-aware anything, unlike scVI; no Rust dependency
or alpha estimation, unlike the archived scclr-based 2.0.concatenate_adata.py). This is meant as
a clean, easy-to-reason-about "before" baseline for comparison against scVI's batch-corrected
results, replacing the old scclr baseline whose OAK outputs were lost (see project memory).

Because this uses plain scanpy PCA (not scclr's all-genes sparse PCA), it uses the standard
HVG-subset -> scale -> PCA -> transfer-embedding-back-onto-full-gene-object pattern.

Each mt threshold (15/20/25/30%) is a fully independent analysis (own normalization, HVG, PCA,
UMAP, Leiden), same as the other step-2 scripts. Final output format matches PerturbNMF's
expected input contract (see 2.1.integrate_scvi.py's docstring for the verified schema).
"""
import gc
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.append("/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb/src")
import utils

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
prepared_scratch_path = os.path.join(os.environ["SCRATCH"], "tf_perturb/2.qc_gex/2.0.prepared_data")
output_path = os.path.join(core_path, "2.qc_gex/2.2.concatenated_data")
scratch_path = os.path.join(os.environ["SCRATCH"], "tf_perturb/2.qc_gex/2.2.concatenated_data")
ensembl_to_symbol_dict = np.load(
    "/oak/stanford/groups/engreitz/Users/opushkar/genome/ensembl_to_symbol_dict_v43.npy",
    allow_pickle=True,
).item()

mt_thresholds = [15, 20, 25, 30]
save_plots = True
random_seed = 0

n_top_genes = 6000
batch_key = "diff_day"      # HVG selection covariate only -- no batch correction happens here
scale_max_value = 10
n_pcs = 50
n_neighbors = 15
leiden_resolutions = [round(r, 2) for r in np.arange(0, 1.3, 0.1)]
leiden_diagnostic_resolution = 0.5
qc_covariates = ["total_counts", "pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"]

labeling_dict = {
    "total_counts": "# UMIs per cell",
    "pct_counts_mt": "% mitochondrial counts",
    "pct_counts_ribo": "% ribosomal counts",
    "S_score": "S score",
    "G2M_score": "G2M score",
    "diff_day": "Differentiation day",
    "rep": "Replicate",
    "well": "Well barcode",
    "guide_category": "Guide category",
    "n_genes_by_counts": "# genes per cell",
}

os.makedirs(scratch_path, exist_ok=True)
for sub in ["plots", "final_data", "tables"]:
    os.makedirs(os.path.join(output_path, sub), exist_ok=True)

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["font.size"] = 14
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.transparent"] = True

np.random.seed(random_seed)

# %%
# One Slurm array task = one mt threshold; falls back to the full sequential loop when run
# directly (no SLURM_ARRAY_TASK_ID).
array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
thresholds_to_run = [mt_thresholds[int(array_task_id)]] if array_task_id is not None else mt_thresholds

for mt_threshold in thresholds_to_run:
    mt_tag = f"mt{mt_threshold}"
    print(f"\n=== {mt_tag} ===")

    checkpoint_path = os.path.join(prepared_scratch_path, f"{mt_tag}_prepared_gex_and_guide.h5mu")
    assert os.path.exists(checkpoint_path), f"Missing 2.0.prepare_adata.py checkpoint: {checkpoint_path}"
    mdata = mu.read(checkpoint_path)
    adata_t, guide_t = mdata["GEX"], mdata["guide"]
    print(f"[{mt_tag}] Loaded {adata_t.n_obs} cells x {adata_t.n_vars} genes")

    # Global well-color palette (all ~96 wells globally unique across the 8 tags)
    all_wells = sorted(adata_t.obs["well"].unique())
    well_palette = np.vstack([plt.cm.tab20.colors, plt.cm.tab20b.colors, plt.cm.tab20c.colors])
    well_colors = {w: well_palette[i % len(well_palette)] for i, w in enumerate(all_wells)}

    # Joint normalization: plain normalize_total (median depth) + log1p, applied once globally
    # -- .layers["counts"] already holds raw counts from 2.0.prepare_adata.py.
    sc.pp.normalize_total(adata_t, target_sum=None)
    sc.pp.log1p(adata_t)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].hist(np.asarray(adata_t.layers["counts"].sum(1)).flatten(), bins=100, color="gray")
    axes[0].set_title("Raw total counts")
    axes[0].set_xlabel(labeling_dict["total_counts"])
    axes[1].hist(np.asarray(adata_t.X.sum(1)).flatten(), bins=100, color="gray")
    axes[1].set_title("Normalized (median-depth + log1p)")
    axes[1].set_xlabel("Normalized total counts")
    plt.tight_layout()
    if save_plots:
        utils.save_fig(os.path.join(output_path, "plots"), f"normalization_check_{mt_tag}")
        pd.DataFrame({
            "raw_total_counts": np.asarray(adata_t.layers["counts"].sum(1)).flatten(),
            "normalized_total_counts": np.asarray(adata_t.X.sum(1)).flatten(),
        }).to_csv(os.path.join(output_path, "plots", f"normalization_check_{mt_tag}_data.tsv"), sep="\t", index=False)
    plt.show()

    # Feature selection: standard "seurat" flavor is appropriate here since .X is genuinely
    # log-normalized (unlike the raw-counts situation in the scclr/scVI scripts).
    sc.pp.highly_variable_genes(adata_t, n_top_genes=n_top_genes, batch_key=batch_key, subset=False)
    ax = sc.pl.highly_variable_genes(adata_t, show=False)
    if save_plots:
        utils.save_fig(os.path.join(output_path, "plots"), f"highly_variable_genes_{mt_tag}")
        adata_t.var[["means", "dispersions", "dispersions_norm", "highly_variable"]].to_csv(
            os.path.join(output_path, "plots", f"highly_variable_genes_{mt_tag}_data.tsv"), sep="\t"
        )
    plt.show()

    adata_hvg = adata_t[:, adata_t.var["highly_variable"]].copy()
    print(f"[{mt_tag}] Subset to {adata_hvg.n_vars} highly variable genes (PCA/clustering input only "
          f"-- the final deliverable below keeps all genes)")

    sc.pp.scale(adata_hvg, max_value=scale_max_value)
    sc.pp.pca(adata_hvg, svd_solver="arpack", random_state=random_seed)
    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=n_pcs, log=True, show=False)
    if save_plots:
        utils.save_fig(os.path.join(output_path, "plots"), f"pca_variance_ratio_{mt_tag}")
        pd.DataFrame({"variance_ratio": adata_hvg.uns["pca"]["variance_ratio"]}).to_csv(
            os.path.join(output_path, "plots", f"pca_variance_ratio_{mt_tag}_data.tsv"), sep="\t", index_label="pc"
        )
    plt.show()

    sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
    sc.tl.umap(adata_hvg, random_state=random_seed)

    # Transfer embedding/neighbors from the HVG subset back onto the full-gene object
    adata_t.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata_t.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata_t.uns["neighbors"] = adata_hvg.uns["neighbors"]
    adata_t.obsp["distances"] = adata_hvg.obsp["distances"]
    adata_t.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
    del adata_hvg

    for res in leiden_resolutions:
        sc.tl.leiden(
            adata_t,
            key_added=f"leiden_res_{res}",
            resolution=res,
            random_state=random_seed,
            flavor="igraph",
            n_iterations=2,
        )

    # UMAP QC panel: day, rep, QC covariates, well, non-targeting guides, side by side
    diag_key = f"leiden_res_{leiden_diagnostic_resolution}"
    panel_keys = ["diff_day", "rep"] + qc_covariates + ["well", "guide_category", diag_key]
    n_cols = 3
    n_rows = int(np.ceil(len(panel_keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.5 * n_cols, 5 * n_rows))
    axes = np.atleast_1d(axes).flatten()
    for i, key in enumerate(panel_keys):
        title = labeling_dict.get(key, key) if key != diag_key else f"Leiden (res={leiden_diagnostic_resolution})"
        palette = well_colors if key == "well" else None
        legend_loc = "none" if key == "well" else "right margin"
        sc.pl.umap(
            adata_t, color=key, ax=axes[i], show=False, title=title,
            palette=palette, legend_loc=legend_loc,
        )
    for ax in axes[len(panel_keys):]:
        ax.axis("off")
    plt.suptitle(f"Simple concat + joint normalization, mt threshold: {mt_threshold}%, n={adata_t.n_obs}")
    plt.tight_layout()
    if save_plots:
        utils.save_fig(os.path.join(output_path, "plots"), f"umap_qc_and_replicate_check_{mt_tag}")
        adata_t.obs[panel_keys + [f"leiden_res_{r}" for r in leiden_resolutions]].to_csv(
            os.path.join(output_path, "plots", f"umap_qc_and_replicate_check_{mt_tag}_data.tsv"), sep="\t"
        )
    plt.show()

    # Replicate AND well composition checks -- both tracked given the ongoing scVI covariate
    # investigation (see project memory).
    for bk in ["rep", "well"]:
        comp_df, overall_chi2, overall_p, dof = utils.replicate_composition_test(
            adata_t, cluster_key=diag_key, batch_key=bk
        )
        comp_df.to_csv(
            os.path.join(output_path, "tables", f"{bk}_composition_test_{mt_tag}.tsv"), sep="\t", index=False
        )
        n_flagged = (comp_df["p_fdr_bh"] < 0.05).sum()
        print(f"[{mt_tag}] {bk} x cluster chi2={overall_chi2:.1f}, dof={dof}, "
              f"FDR-significant: {n_flagged}/{comp_df.shape[0]}")

    # Gene symbol mapping
    adata_t.var["ensembl_id"] = adata_t.var.index.copy()
    adata_t.var["gene_symbol"] = (
        adata_t.var["ensembl_id"].str.split(".").str[0].map(ensembl_to_symbol_dict)
    )
    adata_t.var["gene_symbol"] = adata_t.var["gene_symbol"].fillna(
        adata_t.var["ensembl_id"]
    ).astype(str)
    adata_t.var_names = adata_t.var["gene_symbol"].values
    adata_t.var_names_make_unique()
    adata_t.var.index.name = None

    checkpoint_out = os.path.join(scratch_path, f"concat_normalize_{mt_tag}_all_genes.h5ad")
    adata_t.write_h5ad(checkpoint_out)
    print(f"[{mt_tag}] Wrote scratch checkpoint: {checkpoint_out}")

    # Final deliverable: single AnnData matching PerturbNMF's input contract (see
    # 2.1.integrate_scvi.py's docstring for the verified schema) -- all genes, guide info
    # folded into obsm/uns instead of a separate MuData modality, to $OAK.
    final_barcodes = adata_t.obs_names.intersection(guide_t.obs_names)
    adata_final = adata_t[final_barcodes].copy()
    guide_final = guide_t[final_barcodes].copy()
    adata_final.obs["batch"] = adata_final.obs["diff_day"].values
    adata_final.obsm["guide_assignment"] = sp.csr_matrix(guide_final.layers["guide_assignment"])
    adata_final.uns["guide_names"] = np.array(guide_final.var_names)
    adata_final.uns["guide_targets"] = np.array(guide_final.var["intended_target_name"])
    adata_final.write_h5ad(
        os.path.join(output_path, "final_data", f"concat_normalize_{mt_tag}_all_genes.h5ad")
    )
    print(f"[{mt_tag}] Written final deliverable: {len(final_barcodes)} cells x {adata_final.n_vars} genes")

    del adata_t, adata_final, guide_t, guide_final
    gc.collect()
# %%
