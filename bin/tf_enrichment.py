import pandas as pd
import numpy as np
import pyranges as pr
from pybiomart import Dataset
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from typing import Literal


def run_tf_enrichment(
    results_df: pd.DataFrame,
    peak_file_paths_df: dict[str, str],
    genes: list[str] | None = None,
    element_col: str = "element",
    pval_col: str = "p_value",
    method_col: str = "method",
    gene_col: str = "gene",
    use_gene_ids: bool | Literal["auto"] = "auto",
    fdr_corrected: bool = False,
    promoter_window_width: int = 5000,
    chipseq_threshold: float = 0.75,
    fdr_cutoff: float = 0.1,
    reference: Literal["hg19", "hg38"] = "hg38",
):
    """
    Run TF target enrichment analysis.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame with method results. Required columns: gene_col, method_col, element_col, pval_col.
    peak_file_paths_df : dict
        Dictionary mapping TF names to their ChIP-seq peak file paths.
    genes : list, optional
        List of gene names to consider as potential targets. If None, all genes are used.
    element_col : str, optional
        Column name in results_df that contains TF names. Default is "element".
    pval_col : str, optional
        Column name in results_df that contains p-values. Default is "p_value".
    fdr_corrected : bool, optional
        Whether the p-values in results_df are already FDR corrected. Default is False.
    fdr_cutoff : float, optional
        FDR cutoff for significance. Default is 0.1.
    """
    print(f"Starting TF enrichment analysis using reference genome {reference}...")

    if use_gene_ids == "auto":
        use_gene_ids = any(results_df[gene_col].str.startswith("ENSG"))

    tfs = list(peak_file_paths_df.keys())
    for tf in tfs:
        if tf not in results_df[element_col].values:
            raise ValueError(f"TF {tf} not found in results DataFrame.")

    chipseq_gr = process_chipseq_data(
        peak_file_paths_df,
        element_col=element_col,
        threshold=chipseq_threshold,
    )

    for tf in tfs:
        if chipseq_gr[chipseq_gr.as_df()[element_col] == tf].empty:
            raise ValueError(f"No peaks found for TF {tf} in ChIP-seq data.")

    promoter_gr = get_tss_info(
        gene_names=genes,
        use_ids=use_gene_ids,
        promoter_window_width=promoter_window_width,
        gene_col=gene_col,
        reference=reference,
    )
    # print(promoter_gr)
    # print(chipseq_gr.as_df().head())
    target_matrix = find_tf_targets(
        chipseq_gr, promoter_gr, tfs, element_col=element_col, gene_col=gene_col
    )
    tf_targets = target_matrix.reset_index(names=gene_col).melt(
        id_vars=gene_col, var_name=element_col, value_name="target"
    )

    df_final_input = results_df[results_df[element_col].isin(tfs)].merge(
        tf_targets, how="inner", on=[element_col, gene_col]
    )
    final_summarized_df = summarize_results(
        df_final_input,
        element_col=element_col,
        method_col=method_col,
        fdr_corrected=fdr_corrected,
        pval_col=pval_col,
        q=fdr_cutoff,
    )
    return final_summarized_df


# def compute_gene_by_element_counts(mdata):
#     """Compute gene x element counts."""
#     return pd.DataFrame(
#         mdata["rna"].X.T @ mdata["grna"].X @ mdata["grna"].varm["guide_targets"].values,
#         index=mdata["rna"].var_names,
#         columns=mdata["grna"].varm["guide_targets"].columns,
#     )


# def filter_test_pairs(gene_by_element_counts, mdata, threshold=7, gene_col="gene"):
#     """Filter test pairs based on effective sample size threshold."""
#     control_counts = (
#         mdata["rna"].X.T.astype(bool)
#         @ mdata["grna"][:, mdata["grna"].var_names.str.startswith("NT")].X
#     ).sum(axis=1)
#     passing_genes = mdata["rna"].var_names[control_counts >= threshold]
#     filtered_test_pairs = (
#         gene_by_element_counts.loc[passing_genes, :]
#         .melt(var_name="element", value_name="count", ignore_index=False)
#         .reset_index()
#         .rename(columns={"gene_symbol": gene_col})
#         .query("count >= @threshold")
#     )
#     return filtered_test_pairs


def process_chipseq_data(file_paths, element_col="TF", threshold=0.75):
    """Process ChIP-seq data for transcription factors."""
    col_names = [
        "Chromosome",
        "Start",
        "End",
        "pval",
        "score",
        "pos_max_peak",
        "max_peak_height",
        "rel_pos_max_peak_height",
        "peak_size",
        "mid_point",
        "peak_to_mid_dist",
    ]
    if threshold is None:
        threshold = 0.0
    dfs = []
    for tf, file_path in file_paths.items():
        if file_path.endswith("bed.gz"):
            score_col = "Score"
            df = pr.read_bed(file_path).as_df()
            assert score_col in df.columns, (
                f"Score column {score_col} not found in BED file for TF {tf} among {df.columns}"
            )
        else:
            score_col = "score"
            df = pd.read_csv(file_path, sep="\t", names=col_names, comment="#")
        print(f"Processing TF {tf} with {len(df)} peaks.")
        # df["pval"] = df["pval"].fillna(1)
        # print(df[score_col].quantile(threshold))
        df = df[df[score_col] >= df[score_col].quantile(threshold)]
        print(f"Retained {len(df)} peaks after thresholding for TF {tf}.")
        df[element_col] = tf
        dfs.append(df)
    chipseq_df = pd.concat(dfs, ignore_index=True)
    return pr.PyRanges(chipseq_df)


def get_tss_info(
    gene_names=None,
    use_ids=False,
    promoter_window_width=5000,
    gene_col="gene",
    reference="hg19",
):
    """Query Ensembl for TSS positions and define promoter windows."""

    if reference == "hg19":
        dataset = Dataset(
            name="hsapiens_gene_ensembl", host="http://grch37.ensembl.org"
        )
    elif reference == "hg38":
        dataset = Dataset(name="hsapiens_gene_ensembl", host="http://ensembl.org")
    else:
        raise ValueError(f"Unsupported reference genome: {reference}")
    attrs = [
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
    ]
    if use_ids:
        ensembl_gene_key = "ensembl_gene_id"
    else:
        ensembl_gene_key = "external_gene_name"

    tss_info = dataset.query(attributes=attrs, use_attr_names=True)
    if gene_names is not None:
        tss_info = tss_info[tss_info[ensembl_gene_key].isin(gene_names)]

    valid_chr = [str(i) for i in range(1, 23)] + ["X", "Y"]
    tss_info = tss_info[tss_info["chromosome_name"].isin(valid_chr)]
    tss_info["TSS"] = np.where(
        tss_info["strand"] == 1, tss_info["start_position"], tss_info["end_position"]
    )
    tss_info["chr"] = "chr" + tss_info["chromosome_name"]
    tss_df = tss_info[[ensembl_gene_key, "chr", "TSS"]].rename(
        columns={ensembl_gene_key: gene_col}
    )
    keep_genes = tss_df.groupby(gene_col)["chr"].nunique() == 1
    tss_df = tss_df[tss_df[gene_col].isin(keep_genes[keep_genes].index)]
    tss_summary = (
        tss_df.groupby([gene_col, "chr"], as_index=False)["TSS"]
        .mean()
        .round()
        .astype({"TSS": int})
    )
    tss_summary["Start"] = tss_summary["TSS"] - promoter_window_width
    tss_summary["End"] = tss_summary["TSS"] + promoter_window_width
    return pr.PyRanges(tss_summary.rename(columns={"chr": "Chromosome"}))


def find_tf_targets(
    chipseq_gr, promoter_gr, tf_list, element_col="TF", gene_col="gene"
):
    """Find transcription factor targets by overlap."""
    target_matrix = pd.DataFrame(
        False, index=promoter_gr.df[gene_col].unique(), columns=tf_list
    )
    assert element_col in chipseq_gr.as_df().columns, (
        f"Element column {element_col} not found in ChIP-seq data."
    )
    for TF in tf_list:
        print(f"Processing TF: {TF}")
        tf_sites = chipseq_gr[chipseq_gr.as_df()[element_col] == TF]
        if tf_sites.empty:
            raise ValueError(f"No targets found for TF {TF}")
        overlaps = promoter_gr.join(tf_sites)
        print(f"Found {len(overlaps)} target genes for TF {TF}.")
        genes = overlaps.df[gene_col].unique()
        target_matrix.loc[genes, TF] = True
    return target_matrix


def summarize_results(
    df,
    pval_col="q_value",
    method_col="method",
    element_col="TF",
    q=0.1,
    fdr_corrected=True,
):
    """Summarize results with multiple testing correction."""
    if not fdr_corrected:
        df["significant"] = df.groupby([element_col, method_col])[pval_col].transform(
            lambda pvals: multipletests(pvals, alpha=q, method="fdr_bh")[0]
        )
    else:
        df["significant"] = df[pval_col] <= q
    rows = []
    for (tf, method), sub in df.groupby([element_col, method_col]):
        significant = sub["significant"]
        target_flags = sub["target"]
        if target_flags.nunique() < 2 or significant.nunique() < 2:
            or_val, pval = np.nan, np.nan
        else:
            or_val, pval = fisher_exact(
                [
                    [sum(target_flags & significant), sum(~target_flags & significant)],
                    [
                        sum(target_flags & ~significant),
                        sum(~target_flags & ~significant),
                    ],
                ],
                alternative="greater",
            )
        num_rej = int(significant.sum())
        power = int((significant & target_flags).sum()) / int(target_flags.sum())
        rows.append(
            {
                "TF": tf,
                "Method": method,
                "num_rejections": num_rej,
                "odds_ratio": or_val,
                "pvalue": pval,
                "recall": power,
            }
        )
    # for tf, grp in tf_targets.groupby("TF"):
    #     rows.append(
    #         {
    #             "TF": tf,
    #             "Method": "ChIP-seq",
    #             "num_rejections": int(grp["target"].sum()),
    #             "odds_ratio": np.nan,
    #             "pvalue": np.nan,
    #             "fdp": np.nan,
    #             "power": np.nan,
    #         }
    #     )
    return pd.DataFrame(rows)


# Example usage:
# gene_by_element_counts = compute_gene_by_element_counts(mdata)
# filtered_test_pairs = filter_test_pairs(gene_by_element_counts, mdata, effective_sample_size_threshold)
# chipseq_gr = process_chipseq_data({"STAT1": stat1_txt, "IRF1": irf1_txt}, CHIPSEQ_THRESH)
# promoter_gr = get_tss_info(mdata["rna"].var_names.tolist(), PROMOTER_WINDOW_WIDTH)
# tf_targets = find_tf_targets(chipseq_gr, promoter_gr, ["STAT1", "IRF1"])
# results = summarize_results(df, tf_targets, q)
