#!/usr/bin/env python
# %%
"""
Step 2.0: loads all 8 QC-filtered per-day-rep h5mu files, merges in sceptre guide assignment,
concatenates, filters by guide-assignment status, then writes one mt-filtered (15/20/25/30%)
GEX+guide h5mu checkpoint per threshold for 2.1.integrate_scvi.py to consume.

Split out of 2.0.integrate_scvi.py (2026-09-10) so the (cheap relative to scVI training, but
not cheap in absolute terms -- loads all 8 raw tags before any mt filtering) load/concat/
guide-filter work runs once here instead of being repeated in each of 2.1's 4 parallel scVI
array tasks (that repetition is what caused those tasks to OOM at 256GB -- each was rebuilding
the full unfiltered pool from scratch before ever reaching its own mt subset).

No normalization, PCA, HVG selection, or clustering happens here -- purely load + concat +
guide-filter + mt-split. .X stays raw counts throughout (also copied into layers["counts"] for
downstream naming-convention consistency with 2.1).
"""
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from muon import MuData

sys.path.append("/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb/src")
import utils

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
input_path = os.path.join(core_path, "2.qc_gex/1.qc_data/qc_filtered_data")
guide_assignment_path = os.path.join(core_path, "2.qc_gex/1.qc_data/guide_assignment_data")
output_path = os.path.join(core_path, "2.qc_gex/2.0.prepared_data")
# raw-counts, all-genes checkpoints are large, reproducible intermediates -> $SCRATCH, not $OAK
# (OAK quota is at ~95%; see storage table in the Sherlock CLAUDE.md).
scratch_path = os.path.join(os.environ["SCRATCH"], "tf_perturb/2.qc_gex/2.0.prepared_data")

days = ["d0", "d1", "d2", "d3"]
reps = [1, 2]
save_plots = True

# Narrowed from the full [15, 20, 25, 30] sweep to just the two lenient thresholds
# (user decision, 2026-09-11): the rep1-vs-rep2 quality analysis in 3.rep_quality/ treats
# pct_counts_mt as an *outcome*, and a strict mt filter preferentially removes the worse
# replicate's cells, which would mask the very difference being tested. mt25 is the primary
# analysis, mt30 the leniency sensitivity check. Re-run from scratch rather than reusing the
# Sep-10 checkpoints, also per that decision.
mt_thresholds = [25, 30]

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
os.makedirs(os.path.join(output_path, "plots"), exist_ok=True)

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["font.size"] = 14
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.transparent"] = True

# %%
# Load and directly concatenate all day x rep QC-filtered files
gex_list, guide_list = [], []
guide_assignment_cache = {}  # one file per day (both reps pooled), loaded once and reused
for day in days:
    for rep in reps:
        tag = f"{day}_rep{rep}"
        h5mu_path = os.path.join(input_path, f"{tag}_qc_filtered_gex_and_guide.h5mu")
        assert os.path.exists(h5mu_path), f"Missing h5mu: {h5mu_path}"
        mdata = mu.read(h5mu_path)
        gex_adata, guide_adata = mdata["GEX"], mdata["guide"]
        assert "ENS" in gex_adata.var_names[0]

        # Guide assignment (1.2.guide_assignment_sceptre.py) runs as its own pass, independent
        # of this transcriptome QC/concat pass (see project memory) -- merge its result back in
        # here rather than duplicating the large GEX matrix in that step's output. Assignment
        # is pooled per differentiation day (both reps together, for more cells per guide in
        # sceptre's mixture-model fit), so one file covers both of a day's tags -- subsetting
        # to this tag's barcodes below picks out just this rep's cells.
        if day not in guide_assignment_cache:
            guide_assignment_h5ad = os.path.join(guide_assignment_path, f"{day}_guide_assignment.h5ad")
            assert os.path.exists(guide_assignment_h5ad), f"Missing guide assignment: {guide_assignment_h5ad}"
            guide_assignment_cache[day] = sc.read_h5ad(guide_assignment_h5ad)
        assigned = guide_assignment_cache[day][guide_adata.obs_names]  # align to this tag's barcodes
        guide_adata.layers["guide_assignment"] = assigned.layers["guide_assignment"]
        guide_adata.var["intended_target_name"] = assigned.var["intended_target_name"].values

        print(f"  {tag}: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")
        gex_list.append(gex_adata)
        guide_list.append(guide_adata)
del guide_assignment_cache

adata_concat = ad.concat(gex_list, join="outer", fill_value=0)
# merge="same": ad.concat drops .var columns by default even when identical across inputs
# (verified) -- intended_target_name is identical for a given guide across all 8 tags (same
# guide library, deterministic mapping), so "same" is the correct merge semantics here.
guide_concat = ad.concat(guide_list, join="outer", fill_value=0, merge="same")
del gex_list, guide_list
adata_concat.obs["diff_day"] = adata_concat.obs["day"].astype(str)
print(f"Pooled: {adata_concat.n_obs} cells x {adata_concat.n_vars} genes")
print(f"  day x rep composition: {adata_concat.obs.groupby(['diff_day', 'rep']).size().to_dict()}")

# %%
# Guide-assignment cell filter (independent of mt threshold, applied once)
common_barcodes = adata_concat.obs_names.intersection(guide_concat.obs_names)
adata_concat = adata_concat[common_barcodes].copy()
guide_concat = guide_concat[common_barcodes].copy()

# Drop cells with zero assigned guides (sceptre assignment from 1.2.guide_assignment_sceptre.py,
# merged in above) -- following the pattern of doing guide assignment as its own pass and only
# using it to filter once both the transcriptome and guide data have been aggregated across
# libraries (see project memory).
n_guides_assigned = np.asarray((guide_concat.layers["guide_assignment"] > 0).sum(axis=1)).flatten()
guide_concat.obs["n_guides_assigned"] = n_guides_assigned
has_guide = n_guides_assigned > 0
n_before_guide_filter = adata_concat.n_obs
adata_concat = adata_concat[has_guide].copy()
guide_concat = guide_concat[has_guide].copy()
print(f"guide assignment filter (>=1 assigned guide): {n_before_guide_filter} → {adata_concat.n_obs}")

# Non-targeting identification: a cell counts as "non-targeting" only if every guide assigned to
# it maps to the "non-targeting" control target; any real target present makes it "targeting".
# Independent of mt threshold, so computed once here rather than inside the per-threshold loop.
target_lookup = dict(zip(guide_concat.var_names, guide_concat.var["intended_target_name"]))
guide_names_arr = np.array(guide_concat.var_names)
assignment_csr = sp.csr_matrix(guide_concat.layers["guide_assignment"])
guide_category = []
for row_start, row_end in zip(assignment_csr.indptr[:-1], assignment_csr.indptr[1:]):
    assigned_targets = {target_lookup[g] for g in guide_names_arr[assignment_csr.indices[row_start:row_end]]}
    guide_category.append("non-targeting" if assigned_targets == {"non-targeting"} else "targeting")
guide_concat.obs["guide_category"] = guide_category
adata_concat.obs["guide_category"] = guide_concat.obs["guide_category"].values
adata_concat.obs["n_guides_assigned"] = guide_concat.obs["n_guides_assigned"].values
print(f"  guide_category counts: {adata_concat.obs['guide_category'].value_counts().to_dict()}")

# Global well-color palette (all ~96 wells are globally unique across the 8 tags, no reuse,
# so this needs a larger palette than the per-tag 12-well tab20 palette in 1.qc_adata.py)
all_wells = sorted(adata_concat.obs["well"].unique())
well_palette = np.vstack([plt.cm.tab20.colors, plt.cm.tab20b.colors, plt.cm.tab20c.colors])
well_colors = {w: well_palette[i % len(well_palette)] for i, w in enumerate(all_wells)}

# %%
# Per mt threshold: filter, add raw-counts layer, write GEX+guide checkpoint
for mt_threshold in mt_thresholds:
    mt_tag = f"mt{mt_threshold}"
    print(f"\n=== {mt_tag} ===")

    mt_pass = (adata_concat.obs["pct_counts_mt"] <= mt_threshold).values
    adata_t = adata_concat[mt_pass].copy()
    guide_t = guide_concat[mt_pass].copy()
    print(f"[{mt_tag}] mt filter (<={mt_threshold}%): {adata_concat.n_obs} → {adata_t.n_obs}")

    if save_plots:
        utils.umi_gene_scatter_categorical(
            adata_t,
            labeling_dict=labeling_dict,
            experiment=f"all_days_reps_concat_filtered_{mt_tag}",
            feature="well",
            category_colors=well_colors,
            save=True,
            output_path=os.path.join(output_path, "plots"),
        )

    # Explicit raw-counts layer -- .X is already raw counts (never normalized in this script;
    # 2.1.integrate_scvi.py's scVI trains directly on raw counts via its own NB/ZINB likelihood).
    adata_t.layers["counts"] = adata_t.X.copy()

    checkpoint_path = os.path.join(
        scratch_path, f"{mt_tag}_prepared_gex_and_guide.h5mu")
    mdata_out = MuData({"GEX": adata_t, "guide": guide_t})
    mdata_out.write(checkpoint_path)
    print(f"[{mt_tag}] Wrote: {checkpoint_path}")
# %%
