#!/usr/bin/env python

import mudata as md
from sklearn.metrics import precision_recall_curve, roc_curve, auc
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os


def savefig(path):
    """Utility to save figures cleanly."""
    plt.savefig(path, dpi=300, bbox_inches="tight")
    print(f"✔ Saved plot:", path)


def perform_binary_evaluation(label_controls, infered_significance_col, outdir, plot=True):
    true_label = label_controls
    pred_value = infered_significance_col

    pre, rec, _ = precision_recall_curve(true_label, 1-pred_value)
    auprc = auc(rec, pre)

    fpr, tpr, _ = roc_curve(true_label, 1-pred_value)
    auroc = auc(fpr, tpr)

    print(f"Area under Precision-Recall Curve : {auprc:.3f}")
    print(f"Area under Receiver-Operating Curve: {auroc:.3f}")

    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=130)
        axes[0].plot(rec, pre, lw=1, label=f"AUPRC={auprc:.3f}")
        axes[1].plot(fpr, tpr, lw=1, label=f"AUROC={auroc:.3f}")

        axes[0].set_xlabel("Recall")
        axes[0].set_ylabel("Precision")
        axes[1].set_xlabel("False Positive Rate")
        axes[1].set_ylabel("True Positive Rate")
        axes[0].legend()
        axes[1].legend()
        plt.tight_layout()

        # Save
        savefig(os.path.join(outdir, "trans_perturbo_precision_recall_roc.png"))
        plt.show()


def plot_volcano(
    table_to_fdr,
    outdir,
    fc_col="log2_fc",
    p_col="p_value",
    direct_col="direct_target",
    figsize=(6, 8),
    fc_threshold=1,
    p_threshold=0.05,
    dpi=180
):
    df = table_to_fdr.copy()
    df["neglog10p"] = -np.log10(df[p_col].replace(0, np.nan))

    df["color"] = df[direct_col].map({1: "#d62728", 0: "#7f7f7f"})

    plt.figure(figsize=figsize, dpi=dpi)
    plt.scatter(
        df[fc_col],
        df["neglog10p"],
        c=df["color"],
        s=28,
        alpha=0.75,
        edgecolors="none"
    )

    plt.axvline(x=fc_threshold, color="black", linestyle="--", lw=1, alpha=0.6)
    plt.axvline(x=-fc_threshold, color="black", linestyle="--", lw=1, alpha=0.6)
    plt.axhline(y=-np.log10(p_threshold), color="black", linestyle="--", lw=1, alpha=0.6)

    plt.xlabel("log2 Fold Change", fontsize=13)
    plt.ylabel("-log10(p-value)", fontsize=13)
    plt.title("Volcano Plot", fontsize=15, weight="bold")

    plt.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.35)
    plt.gca().set_facecolor("#f7f7f7")
    for spine in ["top", "right"]:
        plt.gca().spines[spine].set_visible(False)

    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color="#d62728", label="Guides Direct target")
    gray_patch = mpatches.Patch(color="#7f7f7f", label="Non-targeting")
    plt.legend(handles=[red_patch, gray_patch], frameon=False)

    plt.tight_layout()

    # Save
    savefig(os.path.join(outdir, "trans_perturbo_volcano_plot.png"))
    plt.show()


def run_evaluation_controls(md_read, outdir):
    os.makedirs(outdir, exist_ok=True)

    col_used = 'trans_per_guide_results'
    #converting to avoid non boolean values
    col = md_read['guide'].var['targeting']

    md_read['guide'].var['targeting'] = col.apply(
        lambda x: True if (x is True or str(x).upper() == "TRUE")
        else False if (x is False or str(x).upper() == "FALSE")
        else x
    )
        
    selecting_non_targeting_guides = md_read['guide'].var[md_read['guide'].var['targeting'] == False]
    non_targeting_ids = set(selecting_non_targeting_guides['guide_id'].values)

    intended_dict = md_read['guide'].var.set_index(['guide_id'])['intended_target_name'].to_dict()
    md_read.uns[col_used]['intended_target_name'] = md_read.uns[col_used]['guide_id'].map(intended_dict)

    selecting_to_plot = md_read.uns[col_used].copy()
    selecting_to_plot_final = selecting_to_plot[
        selecting_to_plot.apply(lambda x: x['gene_id'] == x['intended_target_name'], axis=1)
    ].drop_duplicates()

    dict_targeting_or_no = md_read['guide'].var['targeting'].to_dict()
    selecting_to_plot['targeting_genes'] = selecting_to_plot['guide_id'].map(dict_targeting_or_no)

    all_targets = selecting_to_plot.query('targeting_genes == True')['intended_target_name'].drop_duplicates()

    non_target_controls = selecting_to_plot.query("targeting_genes == False")
    non_target_controls = non_target_controls[
        non_target_controls.apply(lambda x: x['gene_id'] in all_targets.values, axis=1)
    ]
    non_target_controls['direct_target'] = 0

    table_to_test_cis = md_read.uns[col_used][
        md_read.uns[col_used].apply(lambda x: x['gene_id'] == x['intended_target_name'], axis=1)
    ].drop_duplicates()
    table_to_test_cis['direct_target'] = 1

    table_to_fdr = pd.concat([
        table_to_test_cis,
        non_target_controls.sample(n=table_to_test_cis.shape[0], random_state=42)
    ])

    # Volcano plot
    plot_volcano(table_to_fdr, outdir=outdir)

    # Binary evaluation curves
    perform_binary_evaluation(
        table_to_fdr['direct_target'],
        table_to_fdr['p_value'],
        outdir=outdir,
        plot=True
    )

    # Bar plot
    plt.figure(figsize=(5, 4), dpi=150)
    table_to_fdr.groupby('direct_target').count()['guide_id'].plot(kind='bar')
    plt.ylabel('Number of guides')
    plt.title("Direct targets vs random control guides")
    plt.tight_layout()

    savefig(os.path.join(outdir, "trans_perturbo_barplot_direct_vs_control.png"))
    plt.show()


if __name__ == "__main__":
    print("running controls evaluation program...")
    parser = argparse.ArgumentParser(description="Running controls evaluation program")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--outdir", type=str, default="plots", help="Directory to save plots")

    args = parser.parse_args()

    print("Loading MuData file...")
    
    md_read = md.read_h5mu(args.mdata_path)

    print("Finished Loading MuData file...")

    run_evaluation_controls(md_read, outdir=args.outdir)
