#!/usr/bin/env python

import argparse
import gc

import mudata as mu
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from analysis_output_formatting import (
    format_element_output,
    format_guide_output,
    make_h5mu_safe_dataframe,
)
from mudata_uns_io import write_parquet_dataframe_to_uns, write_uns_patch
from result_table_io import read_result_table, write_result_table
from fast_result_enrichment import (
    enrich_global_parquet,
    supports_fast_parquet_path,
)


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
    print("Loading local-analysis input files...", flush=True)
    local_per_guide = read_result_table(local_analysis_per_guide_path)
    local_per_element = read_result_table(local_analysis_per_element_path)

    local_per_guide = _finalize_perturbo_columns(_add_perturbo_columns(local_per_guide))
    local_per_element = _finalize_perturbo_columns(_add_perturbo_columns(local_per_element))

    # backed="r" avoids loading gene/guide .X into memory; format_guide_output/
    # format_element_output only need .var and .layers["guide_assignment"],
    # which load eagerly even in backed mode.
    base_mdata = mu.read_h5mu(base_mudata_path, backed="r")
    local_per_guide = format_guide_output(local_per_guide, base_mdata)
    local_per_element = format_element_output(local_per_element, base_mdata)
    obsolete_keys = (
        "per_guide_results",
        "per_element_results",
        "cis_per_guide_results",
        "cis_per_element_results",
        "trans_per_guide_results",
        "trans_per_element_results",
    )

    extension = "parquet" if results_format == "parquet" else "tsv.gz"
    write_result_table(local_per_guide, f"local_analysis_per_guide_output.{extension}")
    write_result_table(local_per_element, f"local_analysis_per_element_output.{extension}")

    fast_path = results_format == "parquet" and supports_fast_parquet_path(
        global_analysis_per_guide_path,
        global_analysis_per_element_path,
    )
    if fast_path:
        global_guide_output = "global_analysis_per_guide_output.parquet"
        global_element_output = "global_analysis_per_element_output.parquet"
        try:
            enrich_global_parquet(
                global_analysis_per_guide_path,
                global_guide_output,
                base_mdata,
                "guide",
            )
            enrich_global_parquet(
                global_analysis_per_element_path,
                global_element_output,
                base_mdata,
                "element",
            )
        except ValueError as exc:
            print(
                f"Fast Parquet enrichment is not compatible with these inputs: {exc}. "
                "Falling back to pandas.",
                flush=True,
            )
            fast_path = False

    if not fast_path:
        print("Loading global-analysis inputs with pandas...", flush=True)
        global_per_guide = read_result_table(global_analysis_per_guide_path)
        global_per_element = read_result_table(global_analysis_per_element_path)
        global_per_guide = _finalize_perturbo_columns(
            _add_perturbo_columns(global_per_guide)
        )
        global_per_element = _finalize_perturbo_columns(
            _add_perturbo_columns(global_per_element)
        )
        global_per_guide = format_guide_output(global_per_guide, base_mdata)
        global_per_element = format_element_output(global_per_element, base_mdata)
        write_result_table(
            global_per_guide, f"global_analysis_per_guide_output.{extension}"
        )
        write_result_table(
            global_per_element, f"global_analysis_per_element_output.{extension}"
        )

    base_mdata.file.close()

    # Copy the assay matrices once, then patch one result table at a time. This
    # preserves the self-contained H5MU contract while avoiding simultaneous
    # in-memory copies of both 100M+-row global tables.
    print(f"Writing merged MuData to {output_path}...", flush=True)
    local_guide_safe = make_h5mu_safe_dataframe(local_per_guide)
    local_element_safe = make_h5mu_safe_dataframe(local_per_element)
    write_uns_patch(
        base_mudata_path,
        output_path,
        updates={
            "local_analysis_per_guide_results": local_guide_safe,
            "local_analysis_per_element_results": local_element_safe,
        },
        deletes=obsolete_keys,
    )

    local_guide_rows = len(local_per_guide)
    local_element_rows = len(local_per_element)
    del local_guide_safe, local_element_safe, local_per_guide, local_per_element
    gc.collect()

    if fast_path:
        import pyarrow.parquet as pq

        global_guide_rows = pq.ParquetFile(
            "global_analysis_per_guide_output.parquet"
        ).metadata.num_rows
        write_parquet_dataframe_to_uns(
            output_path,
            "global_analysis_per_guide_results",
            "global_analysis_per_guide_output.parquet",
        )
        gc.collect()

        global_element_rows = pq.ParquetFile(
            "global_analysis_per_element_output.parquet"
        ).metadata.num_rows
        write_parquet_dataframe_to_uns(
            output_path,
            "global_analysis_per_element_results",
            "global_analysis_per_element_output.parquet",
        )
        gc.collect()
    else:
        global_guide_rows = len(global_per_guide)
        global_element_rows = len(global_per_element)
        write_uns_patch(
            output_path,
            output_path,
            updates={
                "global_analysis_per_guide_results": make_h5mu_safe_dataframe(
                    global_per_guide
                )
            },
        )
        del global_per_guide
        gc.collect()
        write_uns_patch(
            output_path,
            output_path,
            updates={
                "global_analysis_per_element_results": make_h5mu_safe_dataframe(
                    global_per_element
                )
            },
        )
        del global_per_element
        gc.collect()

    print("Successfully merged local and global analysis results.", flush=True)
    print(f"  - local_analysis_per_guide_results: {local_guide_rows} entries")
    print(f"  - local_analysis_per_element_results: {local_element_rows} entries")
    print(f"  - global_analysis_per_guide_results: {global_guide_rows} entries")
    print(f"  - global_analysis_per_element_results: {global_element_rows} entries")


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
