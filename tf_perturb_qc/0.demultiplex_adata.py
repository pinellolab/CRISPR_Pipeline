#!/usr/bin/env python
# %%
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from muon import MuData

def round_counts(adata):
    """Round .X to the nearest integer, in place, without densifying.

    kb-count's EM step resolves multi-mapping reads/UMIs probabilistically across
    equivalence classes, so counts_unfiltered/adata.h5ad is not integer-valued (verified
    2026-09-07: e.g. min nonzero GEX value 0.00174216). Several downstream steps model counts
    as draws from a discrete distribution and document/assume true UMI counts (sceptre's
    import_data(): "a matrix of response/gRNA UMI counts"; scVI's default NB likelihood) --
    feeding them fractional values silently degrades the fit rather than erroring. QC
    thresholding, PFlog1pPF, and Scrublet (synthetic_doublet_umi_subsampling=1.0, pure
    addition) don't care either way, so rounding once here, upstream of everything, is
    simpler and safer than remembering to round at each count-distribution-based consumer.
    """
    adata.X = adata.X.tocsr()
    adata.X.data = np.round(adata.X.data)
    adata.X.eliminate_zeros()

sys.path.append("/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb/src")
import utils

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
scratch_path = "/scratch/users/opushkar/full_tf_perturbseq_demux"
sample_info_path = os.path.join(
    core_path, "0.data_preparation/1.input_data/full_TF_sample_info_v1.xlsx"
)
path_to_cc_genes = "/oak/stanford/groups/engreitz/Users/opushkar/common_sc"
parsebio_barcodes_path = os.path.join(
    core_path, "0.data_preparation/1.input_data/bc_data_n198_v5_fixed.csv"
)
output_path = os.path.join(core_path, "2.qc_gex/0.demux_data")
# No cell-calling filter here anymore -- this step just pools raw barcodes across all 13
# sub-libraries per day/rep tag. Cell calling (BarcodeRanks Inflection) happens exactly once,
# at tag scale, in 1.qc_adata.py. Previously this script also applied a per-sub
# BarcodeRanks_Knee filter (from the upstream perturb_pipeline's per-sub cell-calling output)
# before pooling -- removed after verifying (2026-09-04, see project memory) that computing
# BarcodeRanks fresh on the true raw pooled barcodes (1,437,696/tag = 13 x 110,592) gives
# numerically identical thresholds/cell counts to running it on the old per-sub-Knee-filtered
# pool, for both a high-depth tag (d0_rep1) and the lowest-depth tag (d3_rep1, the riskiest
# case for this check). So the per-sub pre-filter was redundant, not wrong -- but keeping it
# meant cell-calling was implemented twice with inconsistent, confusingly-named intermediate
# outputs (see the deleted `n_inflection_barcodes` column, which actually held Knee counts).

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
os.makedirs(os.path.join(output_path, "plots"), exist_ok=True)
os.makedirs(os.path.join(output_path, "data"), exist_ok=True)

# %%
sample_info = pd.read_excel(sample_info_path)
gex_info = sample_info[sample_info.sample_type == "gex"].copy()
gex_info["sub"] = gex_info["sample"].str.extract(r"sub(\d+)").astype(int)
gex_info["day"] = gex_info["sample"].str.extract(r"group_(d\d)")
gex_info["rep"] = gex_info["sample"].str.extract(r"rep(\d+)").astype(int)

print(f"Loaded {len(gex_info)} GEX samples: {sorted(gex_info['day'].unique())} x {sorted(gex_info['rep'].unique())} x {len(gex_info['sub'].unique())} subs")

bc_data = pd.read_csv(parsebio_barcodes_path, comment="#")
seq_to_bci = dict(zip(bc_data["sequence"], bc_data["bci"]))
seq_to_well = dict(zip(bc_data["sequence"], bc_data["well"]))
seq_to_stype = dict(zip(bc_data["sequence"], bc_data["stype"]))
print(f"Loaded {len(seq_to_bci)} barcode→bci/well/stype mappings (stypes: {sorted(bc_data.stype.unique())})")

# %%
summary_rows = []

for (day, rep), group_df in gex_info.groupby(["day", "rep"]):
    tag = f"{day}_rep{rep}"
    print(f"\nProcessing {tag}: {len(group_df)} subpools")

    gex_pools = {}
    guide_pools = {}

    for _, row in group_df.sort_values("sub").iterrows():
        sample_id = row["sample_id"]
        guide_sample_id = row["paired_guide_sample_id"]
        sub = row["sub"]

        gex_h5ad = os.path.join(scratch_path, sample_id, "kb_all_main_raw", "counts_unfiltered", "adata.h5ad")
        guide_h5ad = os.path.join(scratch_path, guide_sample_id, "kb_guide_main_raw", "counts_unfiltered_modified", "adata.h5ad")

        assert os.path.exists(gex_h5ad), f"Missing GEX h5ad: {gex_h5ad}"
        assert os.path.exists(guide_h5ad), f"Missing guide h5ad: {guide_h5ad}"

        gex_adata = sc.read(gex_h5ad)
        guide_adata = sc.read(guide_h5ad)
        round_counts(gex_adata)
        round_counts(guide_adata)

        # No cell-calling here -- just align the two modalities' barcode universes (kb count's
        # GEX and guide-capture whitelists aren't identical). Cell calling happens once, at
        # tag scale, in 1.qc_adata.py.
        common_barcodes = list(set(gex_adata.obs_names).intersection(guide_adata.obs_names))

        print(f"  sub{sub}: {len(gex_adata)} raw barcodes, {len(common_barcodes)} with GEX/guide match")

        gex_sub = gex_adata[common_barcodes].copy()
        guide_sub = guide_adata[common_barcodes].copy()

        # Decode Round1 barcode identity from the last 8 chars of the cell barcode.
        # bci/well/stype come from the fixed bc_data table, which corrects the R-type
        # (random hexamer) well assignments that were colliding with T-type wells.
        round1_seqs = [bc[-8:] for bc in gex_sub.obs_names]
        n_unmatched = sum(seq not in seq_to_bci for seq in round1_seqs)
        if n_unmatched > 0:
            raise RuntimeError(
                f"{n_unmatched}/{len(round1_seqs)} cells in {sample_id} have a Round1 "
                "barcode sequence not found in the fixed bc_data table"
            )
        bci = [seq_to_bci[seq] for seq in round1_seqs]
        well = [seq_to_well[seq] for seq in round1_seqs]
        stype = [seq_to_stype[seq] for seq in round1_seqs]
        gex_sub.obs["bci"] = bci
        guide_sub.obs["bci"] = bci
        gex_sub.obs["well"] = well
        guide_sub.obs["well"] = well
        gex_sub.obs["stype"] = stype
        guide_sub.obs["stype"] = stype

        gex_sub.obs["sub"] = str(sub)
        guide_sub.obs["sub"] = str(sub)
        gex_sub.obs["sample_id"] = sample_id
        guide_sub.obs["sample_id"] = guide_sample_id
        gex_sub.obs["day"] = day
        gex_sub.obs["rep"] = str(rep)
        guide_sub.obs["day"] = day
        guide_sub.obs["rep"] = str(rep)

        gex_sub.obs_names = [f"{bc}_sub{sub}" for bc in gex_sub.obs_names]
        guide_sub.obs_names = [f"{bc}_sub{sub}" for bc in guide_sub.obs_names]

        gex_pools[sub] = gex_sub
        guide_pools[sub] = guide_sub

        stype_counts = pd.Series(stype).value_counts().to_dict()
        summary_rows.append({
            "tag": tag, "day": day, "rep": rep, "sub": sub,
            "n_raw_barcodes": len(gex_adata),
            "n_common_barcodes": len(common_barcodes),
            "stype_counts": stype_counts,
        })

    gex_concat = ad.concat(list(gex_pools.values()), merge="same")
    guide_concat = ad.concat(list(guide_pools.values()), merge="same")

    print(f"  {tag} concat: {gex_concat.n_obs} raw (GEX/guide-matched) barcodes x {gex_concat.n_vars} genes")

    mdata = MuData({"GEX": gex_concat, "guide": guide_concat})
    out_h5mu = os.path.join(output_path, "data", f"{tag}_all_subs_gex_and_guide.h5mu")
    mdata.write(out_h5mu)
    print(f"  Written: {out_h5mu}")

# %%
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(output_path, "data", "demux_cell_counts_summary.tsv"), sep="\t", index=False)
print("\nDemux summary:")
print(summary_df.groupby(["day", "rep"])[["n_raw_barcodes", "n_common_barcodes"]].sum().to_string())
# %%
