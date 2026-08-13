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
from mudata_uns_io import write_uns_patch
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
        if "perturbo_q_value" not in out.columns:
            out["perturbo_q_value"] = _bh_adjust(out["perturbo_p_value"])
        else:
            missing_q = out["perturbo_q_value"].isna() & out["perturbo_p_value"].notna()
            if missing_q.any():
                computed_q = _bh_adjust(out["perturbo_p_value"])
                out.loc[missing_q, "perturbo_q_value"] = computed_q.loc[missing_q]
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

    # backed="r" avoids loading gene/guide .X into memory; format_guide_output/
    # format_element_output only need .var and .layers["guide_assignment"],
    # which load eagerly even in backed mode.
    base_mdata = mu.read_h5mu(base_mudata_path, backed="r")
    local_per_guide = format_guide_output(local_per_guide, base_mdata)
    global_per_guide = format_guide_output(global_per_guide, base_mdata)
    local_per_element = format_element_output(local_per_element, base_mdata)
    global_per_element = format_element_output(global_per_element, base_mdata)

    print(f"Writing merged MuData to {output_path}...")
    # make_h5mu_safe_dataframe encodes the low-cardinality string columns in
    # these tables (gene/guide ids, target names, chromosomes, pair types) as
    # categoricals, which captures most of the disk-size win gzip would
    # otherwise be relied on for. That makes the fast write_uns_patch path
    # (byte-copy the matrices, patch only /uns) the better trade here instead
    # of a full compressed re-serialize.
    updates = {
        "local_analysis_per_guide_results": make_h5mu_safe_dataframe(local_per_guide),
        "local_analysis_per_element_results": make_h5mu_safe_dataframe(local_per_element),
        "global_analysis_per_guide_results": make_h5mu_safe_dataframe(global_per_guide),
        "global_analysis_per_element_results": make_h5mu_safe_dataframe(global_per_element),
    }
    obsolete_keys = (
        "per_guide_results",
        "per_element_results",
        "cis_per_guide_results",
        "cis_per_element_results",
        "trans_per_guide_results",
        "trans_per_element_results",
    )
    base_mdata.file.close()
    write_uns_patch(
        base_mudata_path, output_path, updates=updates, deletes=obsolete_keys
    )

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
