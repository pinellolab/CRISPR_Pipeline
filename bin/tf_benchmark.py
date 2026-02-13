#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import muon as mu

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from tf_enrichment import run_tf_enrichment


PEAK_META_DICT = {
    "ENSG00000166949": "ENCFF580MYL.bed.gz",  # SMAD3
    "ENSG00000081059": "ENCFF032CCV.bed.gz",  # TCF7
    "ENSG00000141646": "ENCFF126SGU.bed.gz",  # SMAD4
    "ENSG00000149948": "ENCFF614IHQ.bed.gz",  # HMGA2
    "ENSG00000072364": "ENCFF835WQS.bed.gz",  # AFF4
}

ENSEMBL_GENE_NAME_MAP = {
    "ENSG00000166949": "SMAD3",
    "ENSG00000081059": "TCF7",
    "ENSG00000141646": "SMAD4",
    "ENSG00000149948": "HMGA2",
    "ENSG00000072364": "AFF4",
}


def _load_results(mdata):
    results = mdata.uns.get("trans_per_element_results")
    if results is None:
        raise ValueError("Missing trans_per_element_results in MuData.uns")

    if not isinstance(results, pd.DataFrame):
        results = pd.DataFrame(results)

    results = results.drop_duplicates().copy()
    if "intended_target_name" not in results.columns:
        if "guide_id" not in results.columns:
            raise ValueError("Missing guide_id in trans_per_element_results.")
        if "guide" not in mdata.mod:
            raise ValueError("Missing guide modality in MuData.")
        guide_to_target = mdata["guide"].var.get("intended_target_name")
        if guide_to_target is None:
            raise ValueError("Missing intended_target_name in guide var metadata.")
        results = results.assign(
            method="PerTurbo",
            intended_target_name=lambda x: (
                x["guide_id"]
                .map(guide_to_target)
                .astype(str)
                .replace("nan", "non-targeting")
            ),
        )
        results = (
            results.groupby(["gene_id", "intended_target_name", "method"], as_index=False)
            .agg(
                log2_fc=("log2_fc", "mean"),
                p_value=("p_value", "min"),
            )
        )
    else:
        if "method" not in results.columns:
            results = results.assign(method="PerTurbo")
    return results


def _get_gene_universe(mdata):
    if "gene" in mdata.mod:
        return mdata["gene"].var_names.tolist()
    if "rna" in mdata.mod:
        return mdata["rna"].var_names.tolist()
    return None


def _resolve_peak_files(encode_bed_dir):
    paths = {}
    for gene_id, filename in PEAK_META_DICT.items():
        path = os.path.join(encode_bed_dir, filename)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing ENCODE BED file: {path}")
        paths[gene_id] = path
    return paths


def generate_tf_benchmark(
    mudata_path,
    encode_bed_dir,
    output_dir,
    reference_gtf=None,
    dataset_name="Current run",
):
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "benchmark_tables")
    os.makedirs(tables_dir, exist_ok=True)

    mdata = mu.read(mudata_path)
    results = _load_results(mdata)
    genes_universe = _get_gene_universe(mdata)

    peak_file_paths = _resolve_peak_files(encode_bed_dir)
    mapping_df = pd.DataFrame(
        [
            {"tf_id": k, "tf_name": ENSEMBL_GENE_NAME_MAP.get(k, k), "bed_file": v}
            for k, v in peak_file_paths.items()
        ]
    )
    mapping_df.to_csv(os.path.join(tables_dir, "tf_peak_mapping.tsv"), sep="\t", index=False)
    results.to_csv(os.path.join(tables_dir, "trans_per_element_results_used.tsv"), sep="\t", index=False)

    promoter_windows = [500, 1000, 2500, 5000, 10000]

    enrichment_runs = []
    for w in promoter_windows:
        er = run_tf_enrichment(
            results.assign(method="PerTurbo"),
            peak_file_paths,
            gene_col="gene_id",
            element_col="intended_target_name",
            genes=genes_universe,
            chipseq_threshold=0.00,
            promoter_window_width=w,
            fdr_corrected=True,
            fdr_cutoff=0.001,
            reference=reference_gtf if reference_gtf else "hg38",
        )
        if er is None or len(er) == 0:
            continue
        er = er.copy()
        er["TF_display"] = er["TF"].map(ENSEMBL_GENE_NAME_MAP).fillna(er["TF"])
        er["promoter_window_width"] = w
        enrichment_runs.append(er)

    if not enrichment_runs:
        raise ValueError("No enrichment results were produced. Check inputs/paths.")

    enrichment_all = pd.concat(enrichment_runs, ignore_index=True)
    enrichment_all["log_pvalue_clipped"] = np.log10(
        enrichment_all["pvalue"].astype(float).clip(lower=1e-100, upper=1.0)
    )

    enrichment_all.to_csv(
        os.path.join(tables_dir, "enrichment_all.tsv"), sep="\t", index=False
    )

    tf_order = (
        enrichment_all.groupby("TF_display")["pvalue"]
        .min()
        .sort_values(ascending=True)
        .index.tolist()
    )
    pd.DataFrame({"TF": tf_order}).to_csv(
        os.path.join(tables_dir, "tf_order.tsv"), sep="\t", index=False
    )

    sns.set(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    sns.barplot(
        data=enrichment_all,
        x="num_rejections",
        y="TF_display",
        hue="promoter_window_width",
        order=tf_order,
        ax=axes[0],
    )
    axes[0].set_xscale("log")
    axes[0].set_xlim(1, 1e4)
    axes[0].set_xlabel("N downstream genes")
    axes[0].set_ylabel("Transcription factor")

    sns.barplot(
        data=enrichment_all.dropna(subset=["odds_ratio"]),
        x="odds_ratio",
        y="TF_display",
        hue="promoter_window_width",
        order=tf_order,
        ax=axes[1],
    )
    axes[1].axvline(1, linestyle="--", linewidth=1)
    axes[1].set_xlabel("Odds ratio")
    axes[1].set_ylabel("")

    sns.barplot(
        data=enrichment_all.dropna(subset=["log_pvalue_clipped"]),
        x="log_pvalue_clipped",
        y="TF_display",
        hue="promoter_window_width",
        order=tf_order,
        ax=axes[2],
    )
    axes[2].axvline(np.log10(0.05), linestyle="--", linewidth=1)
    axes[2].set_xlabel(r"$\log_{10}$(p-value)")
    axes[2].set_ylabel("")

    axes[0].legend_.remove()
    axes[1].legend_.remove()
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=len(promoter_windows),
        title="Promoter window (bp)",
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle(dataset_name)
    fig_path = os.path.join(output_dir, "tf_benchmark.png")
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Run TF benchmark plots.")
    parser.add_argument("--mudata", required=True, help="Path to inference_mudata.h5mu")
    parser.add_argument("--encode_bed_dir", required=True, help="Path to ENCODE BED files directory")
    parser.add_argument("--output_dir", default="benchmark_output", help="Output directory")
    parser.add_argument("--gtf", default=None, help="Reference GTF (optional)")
    parser.add_argument("--dataset_name", default="Current run", help="Dataset name for plot title")
    args = parser.parse_args()

    generate_tf_benchmark(
        mudata_path=args.mudata,
        encode_bed_dir=args.encode_bed_dir,
        output_dir=args.output_dir,
        reference_gtf=args.gtf,
        dataset_name=args.dataset_name,
    )


if __name__ == "__main__":
    main()
