#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline


def get_vectors(x, y):
    smooth_spline = UnivariateSpline(x, y, s=len(x))
    second_deriv = smooth_spline.derivative(n=2)(x)

    ten_percent = max(1, round(len(x) * 0.1))
    if len(second_deriv) > 2 * ten_percent:
        mid_second_deriv = second_deriv[ten_percent:-ten_percent]
    else:
        mid_second_deriv = second_deriv

    if np.all(mid_second_deriv >= 0) or np.all(mid_second_deriv <= 0):
        return x, y

    abs_min_pos = int(np.argmin(second_deriv))
    left_vect = second_deriv[:abs_min_pos + 1]
    endpt_1_candidates = np.where(left_vect >= 0)[0]
    if len(endpt_1_candidates) == 0:
        return x, y
    endpt_1_idx = endpt_1_candidates[-1]

    right_vect = second_deriv[abs_min_pos:]
    endpt_2_candidates = np.where(right_vect >= 0)[0]
    if len(endpt_2_candidates) == 0:
        return x, y
    endpt_2_idx = abs_min_pos + endpt_2_candidates[0]

    if endpt_1_idx >= endpt_2_idx:
        return x, y

    return x[endpt_1_idx:endpt_2_idx + 1], y[endpt_1_idx:endpt_2_idx + 1]


def elbow_knee_finder(x, y, mode="basic"):
    if mode == "advanced":
        if len(np.unique(x)) < 4:
            return None
        x, y = get_vectors(x, y)

    if len(x) == 0 or len(y) == 0:
        return None

    x0, y0 = x[0], y[0]
    x1, y1 = x[-1], y[-1]

    if x0 == x1:
        return None

    slope = (y1 - y0) / (x1 - x0)
    intercept = y0 - slope * x0
    distances = np.abs(slope * x - y + intercept) / np.sqrt(slope ** 2 + 1)

    max_idx = int(np.argmax(distances))
    return np.array([x[max_idx], y[max_idx]])


def get_elbow_knee_points(x, y):
    point_1 = elbow_knee_finder(x, y, mode="basic")
    point_2 = None
    if point_1 is not None:
        end_idx = int(round(point_1[0]))
        end_idx = max(1, min(len(x), end_idx))
        point_2 = elbow_knee_finder(x[:end_idx], y[:end_idx], mode="advanced")
    return point_1, point_2


def plot_barcode_rank(knee_df, point_1, point_2, outpath):
    rank = knee_df["rank"].values
    log_counts = knee_df["sum_log"].values

    plt.figure(figsize=(8, 5))
    plt.plot(rank, log_counts, marker="o", linestyle="-", markersize=2)

    if point_1 is not None:
        knee_rank = int(round(point_1[0]))
        plt.axvline(knee_rank, color="red", linestyle="--", alpha=0.8, label="Knee 1")
        plt.scatter([point_1[0]], [point_1[1]], color="red", s=30)

    if point_2 is not None:
        knee_rank_2 = int(round(point_2[0]))
        plt.axvline(knee_rank_2, color="orange", linestyle="--", alpha=0.8, label="Knee 2")
        plt.scatter([point_2[0]], [point_2[1]], color="orange", s=30)

    plt.ylabel("Log1p UMI Counts")
    plt.xlabel("Barcode Rank")
    plt.title("Barcode Rank Plot")
    if point_1 is not None or point_2 is not None:
        plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def main(adata_rna, gname_rna, min_genes, min_cells, pct_mito, reference, barcode_filter):

    gname_rna_path = os.path.join(gname_rna, 'counts_unfiltered/cells_x_genes.genes.names.txt')
    gene_df = pd.read_csv(gname_rna_path, header=None)
    gene_names = gene_df[0].tolist()
    adata_rna = sc.read(adata_rna)

    if len(gene_names) == adata_rna.shape[1]:
        # keep ensembl id as var_names
        # matched gene names as column: gene_id
        # adata_rna.var_names = gene_names
        adata_rna.var["symbol"] = gene_names
        # modify the ensembl id in adata_rna.var_names
        adata_rna.var_names = adata_rna.var_names.str.split('.').str[0]
        adata_rna.var_names_make_unique()
    else: 
        raise ValueError("The number of gene names does not match the number of variables in adata_rna")

    # Save knee plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")

    cell_counts = np.asarray(adata_rna.X.sum(axis=1)).ravel()
    knee_df = pd.DataFrame({
        'sum': cell_counts,
        'barcodes': adata_rna.obs_names.values})
    knee_df = knee_df.sort_values('sum', ascending=False).reset_index(drop=True)
    knee_df['sum_log'] = np.log1p(knee_df['sum'])
    knee_df['rank'] = np.arange(1, len(knee_df) + 1)

    point_1, point_2 = get_elbow_knee_points(
        knee_df["rank"].values, knee_df["sum_log"].values
    )
    plot_barcode_rank(knee_df, point_1, point_2, 'figures/knee_plot_scRNA.png')

    if barcode_filter != "none":
        selected_point = point_1 if barcode_filter == "knee" else point_2
        if selected_point is None:
            print("Barcode filtering skipped: knee point could not be determined.")
        else:
            knee_rank = int(round(selected_point[0]))
            knee_rank = max(1, min(len(knee_df), knee_rank))
            count_threshold = knee_df.loc[knee_rank - 1, "sum"]
            print(f"Barcode filtering using rank {knee_rank} (count >= {count_threshold:.0f})")
            keep_mask = cell_counts >= count_threshold
            adata_rna = adata_rna[keep_mask, :].copy()

    # Add batch number
    adata_rna.obs['batch_number'] = adata_rna.obs['batch'].factorize()[0] + 1

    mt_prefix = "MT-" if reference == "human" else "Mt-"
    adata_rna.var["mt"] = adata_rna.var['symbol'].str.startswith('MT-')
    adata_rna.var["ribo"] = adata_rna.var['symbol'].str.startswith(("RPS", "RPL"))

    # Calculate QC metrics
    adata_rna.X = adata_rna.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata_rna, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

    # Plot violin
    print (adata_rna)
    adata_rna
    sc.pl.violin(
        adata_rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save='plot_scrna.png'
    )

    # Plot scatter
    sc.pl.scatter(
        adata_rna,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        size=0.1,
        save='plot_scrna.png'
    )

    # filter for mic_cells and min_genes
    if barcode_filter == "none":
        sc.pp.filter_cells(adata_rna, min_genes=min_genes)
    else:
        print("Skipping min_genes filter because barcode_filter is enabled.")
    sc.pp.filter_genes(adata_rna, min_cells=10) # just to remove some basic amout of genes, this will be treated in during the mudata concat

    # filter for percent mito
    pct_mito=pct_mito 
    adata_rna = adata_rna[adata_rna.obs['pct_counts_mt'] < pct_mito, :]

    # Save
    adata_rna.write('filtered_anndata.h5ad')

    print("Quality control completed and saved to 'filtered_anndata.h5ad'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC on AnnData.')
    parser.add_argument('adata_rna', type=str, help='Path to the AnnData file.')
    parser.add_argument('gname_rna', type=str, help='Path to the cells x genes txt file.')
    parser.add_argument('--min_genes', type=int, default=100, help='Minimum number of genes per cell.')
    parser.add_argument('--min_cells', type=float, default=0.05, help='Minimum number of cells per gene.')
    parser.add_argument('--pct_mito', type=float, default=0.2, help='Minimum percent of proportion of mitochondrial reads in cells.')
    parser.add_argument('--reference', type=str, required=True, help='Reference species')
    parser.add_argument(
        '--barcode-filter',
        choices=['none', 'knee', 'knee2'],
        default='knee',
        help="Barcode filtering method using knee detection (default: 'knee').",
    )

    args = parser.parse_args()
    main(
        args.adata_rna,
        args.gname_rna,
        args.min_genes,
        args.min_cells,
        args.pct_mito,
        args.reference,
        args.barcode_filter,
    )
