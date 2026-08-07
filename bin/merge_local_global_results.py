#!/usr/bin/env python

import argparse

import mudata as mu
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from analysis_output_formatting import (
    format_element_output,
    format_guide_output,
    make_h5mu_safe_dataframe,
)
from result_table_io import read_result_table, write_result_table


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if not valid.any():
        return out

    clipped = p.loc[valid].clip(lower=0.0, upper=1.0)
    out.loc[clipped.index] = false_discovery_control(
        clipped.to_numpy(dtype=float), method="bh"
    )
    return out


def _add_perturbo_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(
        columns={
            "log2_fc": "perturbo_log2_fc",
            "p_value": "perturbo_p_value",
            "log2_fc_std": "perturbo_fc_se",
            "q_value": "perturbo_q_value",
        }
    )
    out = out.drop(columns=["perturbo_fdr_log10_p_value"], errors="ignore")
    if "perturbo_p_value" in out.columns:
        out["perturbo_q_value"] = _bh_adjust(out["perturbo_p_value"])
    return out


def _finalize_perturbo_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "perturbo_fc_se" not in out.columns and "perturbo_p_value" in out.columns:
        out["perturbo_fc_se"] = np.nan
    return out


def merge_local_global_results(
    local_analysis_per_guide_path,
    local_analysis_per_element_path,
    global_analysis_per_guide_path,
    global_analysis_per_element_path,
    base_mudata_path,
    output_path,
    results_format="parquet",
):
    """Merge local- and global-analysis result tables into the final MuData."""
    print("Loading input files...")
    local_per_guide = read_result_table(local_analysis_per_guide_path)
    local_per_element = read_result_table(local_analysis_per_element_path)
    global_per_guide = read_result_table(global_analysis_per_guide_path)
    global_per_element = read_result_table(global_analysis_per_element_path)

    local_per_guide = _finalize_perturbo_columns(_add_perturbo_columns(local_per_guide))
    local_per_element = _finalize_perturbo_columns(_add_perturbo_columns(local_per_element))
    global_per_guide = _finalize_perturbo_columns(_add_perturbo_columns(global_per_guide))
    global_per_element = _finalize_perturbo_columns(_add_perturbo_columns(global_per_element))

    base_mdata = mu.read_h5mu(base_mudata_path)
    local_per_guide = format_guide_output(local_per_guide, base_mdata)
    global_per_guide = format_guide_output(global_per_guide, base_mdata)
    local_per_element = format_element_output(local_per_element, base_mdata)
    global_per_element = format_element_output(global_per_element, base_mdata)

    base_mdata.uns["local_analysis_per_guide_results"] = make_h5mu_safe_dataframe(local_per_guide)
    base_mdata.uns["local_analysis_per_element_results"] = make_h5mu_safe_dataframe(local_per_element)
    base_mdata.uns["global_analysis_per_guide_results"] = make_h5mu_safe_dataframe(global_per_guide)
    base_mdata.uns["global_analysis_per_element_results"] = make_h5mu_safe_dataframe(global_per_element)

    for obsolete_key in (
        "per_guide_results",
        "per_element_results",
        "cis_per_guide_results",
        "cis_per_element_results",
        "trans_per_guide_results",
        "trans_per_element_results",
    ):
        base_mdata.uns.pop(obsolete_key, None)

    print(f"Writing merged MuData to {output_path}...")
    base_mdata.write(output_path, compression="gzip")

    extension = "parquet" if results_format == "parquet" else "tsv.gz"
    write_result_table(local_per_guide, f"local_analysis_per_guide_output.{extension}")
    write_result_table(local_per_element, f"local_analysis_per_element_output.{extension}")
    write_result_table(global_per_guide, f"global_analysis_per_guide_output.{extension}")
    write_result_table(global_per_element, f"global_analysis_per_element_output.{extension}")

    print("Successfully merged local and global analysis results.")
    print(f"  - local_analysis_per_guide_results: {len(local_per_guide)} entries")
    print(f"  - local_analysis_per_element_results: {len(local_per_element)} entries")
    print(f"  - global_analysis_per_guide_results: {len(global_per_guide)} entries")
    print(f"  - global_analysis_per_element_results: {len(global_per_element)} entries")
    return base_mdata


def main():
    parser = argparse.ArgumentParser(description="Merge local- and global-analysis result tables")
    parser.add_argument("--local_analysis_per_guide", required=True)
    parser.add_argument("--local_analysis_per_element", required=True)
    parser.add_argument("--global_analysis_per_guide", required=True)
    parser.add_argument("--global_analysis_per_element", required=True)
    parser.add_argument("--base_mudata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--results_format",
        choices=["tsv.gz", "parquet"],
        default="parquet",
        help="Format for the four merged result tables.",
    )
    args = parser.parse_args()

    merge_local_global_results(
        args.local_analysis_per_guide,
        args.local_analysis_per_element,
        args.global_analysis_per_guide,
        args.global_analysis_per_element,
        args.base_mudata,
        args.output,
        args.results_format,
    )


if __name__ == "__main__":
    main()
