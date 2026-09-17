#!/usr/bin/env python3
"""Concatenate measurement-set-filtered RNA AnnData and apply global gene QC."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from count_matrix_utils import normalize_sparse_index_dtypes


def validate_feature_index(inputs):
    expected = None
    template_var = None
    for path in inputs:
        current = ad.read_h5ad(path, backed="r")
        index = pd.Index(current.var_names.astype(str))
        if template_var is None:
            template_var = current.var.copy()
        current.file.close()
        if expected is None:
            expected = index
        elif not index.equals(expected):
            raise ValueError(
                "QC-filtered measurement sets do not share the same ordered "
                f"RNA feature index: {path} differs from {inputs[0]}"
            )
    return template_var


def recompute_gene_metrics(adata):
    detected_cells = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    total_counts = np.asarray(adata.X.sum(axis=0)).ravel()
    mean_counts = total_counts / adata.n_obs if adata.n_obs else np.zeros(adata.n_vars)
    adata.var["n_cells_by_counts"] = detected_cells
    adata.var["mean_counts"] = mean_counts
    adata.var["log1p_mean_counts"] = np.log1p(mean_counts)
    adata.var["pct_dropout_by_counts"] = (
        100 * (1 - detected_cells / adata.n_obs) if adata.n_obs else np.nan
    )
    adata.var["total_counts"] = total_counts
    adata.var["log1p_total_counts"] = np.log1p(total_counts)
    return detected_cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output", default="filtered_anndata.h5ad")
    parser.add_argument("--tapseq-mode", action="store_true")
    args = parser.parse_args()
    inputs = sorted(map(Path, args.inputs), key=lambda path: path.name)
    if not inputs:
        parser.error("At least one filtered measurement-set AnnData is required")
    template_var = validate_feature_index(inputs)
    ad.experimental.concat_on_disk(
        inputs,
        join="outer",
        # AnnData 0.11 cannot write pandas Series from merge="first" in
        # concat_on_disk. Feature indices are validated above, so restore the
        # first measurement set's feature metadata explicitly after concat.
        merge=None,
        index_unique=None,
        out_file=Path(args.output),
    )
    combined = ad.read_h5ad(args.output)
    combined.var = template_var.loc[combined.var_names].copy()
    if not combined.obs_names.is_unique:
        raise ValueError("Per-measurement-set QC produced duplicate qualified cell barcodes")
    if "batch" in combined.obs:
        combined.obs["batch_number"] = combined.obs["batch"].factorize()[0] + 1
    combined.X = normalize_sparse_index_dtypes(combined.X)
    minimum_cells = 1 if args.tapseq_mode else 10
    detected_cells = recompute_gene_metrics(combined)
    combined = combined[:, detected_cells >= minimum_cells].copy()
    recompute_gene_metrics(combined)
    combined.write_h5ad(args.output)
    print(f"Concatenated {len(inputs)} QC-filtered measurement sets; retained {combined.n_obs} cells and {combined.n_vars} genes")


if __name__ == "__main__":
    main()
