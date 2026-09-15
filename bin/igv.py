#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import mudata as mu
import os
from gtfparse import read_gtf
from typing import Literal, Tuple, Dict, Optional, Any

# Why this file is column-at-a-time rather than row-at-a-time:
#
# ``global_analysis_per_element_results`` on the TAP-seq chr8 screen is
# 13,140,543 rows, and ``evaluation_plot`` walked it with ``DataFrame.iterrows``
# -- one materialised Series per row -- to decide between two branches and
# append five or eight scalars. That task took 1.19 h (trace realtime
# 4,275,689 ms) to write 78,690 output rows.
#
# Every per-row decision here is a boolean mask instead, and the output columns
# are built only for the rows that survive. Two things the row loop got for free
# have to be held deliberately:
#
#  * Branch choice. ``row[a] == row[b]`` is Python ``==``, so ``None == None`` is
#    True while ``nan == nan`` is False. Object-dtype ``==`` on ndarrays is
#    elementwise ``PyObject_RichCompare`` and agrees; pandas' own ``==`` does
#    not (it forces False whenever either side is null).
#  * Output dtype. The bedpe/bedgraph files are pipeline artifacts, and
#    ``to_csv`` formats from dtype. Gene coordinates arrive as float64 and guide
#    coordinates as nullable Int64, so the same column prints ``85107147.0`` or
#    ``94879769`` depending on which dictionary entry fed it; and the global
#    table's log2_fc is float32, which the row loop nonetheless widened to
#    float64. Both fall out of one fact: ``iterrows`` reads cells from
#    ``DataFrame.values``, whose object cast unboxes every numpy scalar. So the
#    values handed to ``pd.DataFrame`` here are the same Python scalars the row
#    loop appended, and inference is left to do what it did before. Nothing may
#    go through ``Series.map`` or reach the frame as a typed ndarray.


def _row_values(var: pd.DataFrame) -> Tuple[np.ndarray, Dict[str, int]]:
    """A ``.var`` frame as the array ``iterrows`` walks, plus column positions.

    ``DataFrame.iterrows`` iterates ``DataFrame.values``, so the scalar a row
    lookup returns depends on the frame's *common* dtype, not the column's: a
    frame carrying any string column is object, and every float64 cell in it
    comes back as a Python float. Taking the same whole-frame array keeps those
    scalars byte-identical downstream, which per-column ``to_numpy`` would not
    -- and it is the per-row Series construction that was slow, not this.
    """
    return var.to_numpy(), {name: i for i, name in enumerate(var.columns)}


def _isnan_mask(values: np.ndarray) -> np.ndarray:
    """``np.isnan`` over a whole column, matching the per-value guard.

    Float columns go through ``np.isnan`` directly. Integer and boolean columns
    can hold no nan, which is what ``np.isnan`` would have said value by value.
    Anything else -- object, which is what a mixed ``.var`` frame gives -- is
    tested one value at a time with the same call, so a column ``np.isnan``
    cannot accept still raises the same ``TypeError`` it used to.
    """
    if values.dtype.kind in "fc":
        return np.isnan(values)
    if values.dtype.kind in "iub":
        return np.zeros(values.shape, dtype=bool)
    return np.array([bool(np.isnan(value)) for value in values], dtype=bool)


def _contains(values: np.ndarray, mapping: Dict[Any, list]) -> np.ndarray:
    """``[v in mapping for v in values]`` as a boolean array.

    The dict probe is kept rather than swapped for ``isin`` or a factorize so
    that nan and None keys resolve exactly the way ``in`` resolved them for the
    row loop. It costs ~0.8 s per 13.1M-row column.
    """
    return np.fromiter(
        (value in mapping for value in values), dtype=bool, count=len(values)
    )


def _selected(series: pd.Series, mask: np.ndarray) -> list:
    """The masked column as the row loop's list of scalars.

    ``tolist``, not the ndarray, and the difference is visible in the output
    file. The merged frame always carries string columns, so it is object, so
    ``iterrows`` read every cell out of an object array -- and an object cast
    unboxes numpy scalars into Python ones. ``global_analysis_per_element_results``
    stores ``perturbo_log2_fc`` as float32; via Python floats that column infers
    back to float64 and prints ``-0.17034076154232025``, while handing over the
    float32 array prints ``-0.17034076``. ``tolist`` performs the same unboxing,
    and leaves object arrays (extension dtypes, pd.NA) untouched.
    """
    return series.to_numpy()[mask].tolist()


def process_coordinates(mdata) -> Dict[str, list]:
    """Extract coordinate information from MuData object"""
    coordinate_dict = {}

    # Process gene coordinates
    gene_var = mdata.mod["gene"].var
    rows, at = _row_values(gene_var)
    gene_start = rows[:, at["gene_start"]]
    gene_end = rows[:, at["gene_end"]]
    placed = ~(_isnan_mask(gene_start) | _isnan_mask(gene_end))
    for index, chrom, start, end in zip(
        gene_var.index.to_numpy()[placed],
        rows[placed, at["gene_chr"]],
        gene_start[placed],
        gene_end[placed],
    ):
        coordinate_dict[index] = [chrom, start, end]

    # Process guide coordinates. Still a sequential pass: each row is skipped
    # against the dictionary as it stands, so gene entries win over guides and
    # the first guide row for a target wins over later ones.
    rows, at = _row_values(mdata.mod["guide"].var)
    for name, chrom, start, end in zip(
        rows[:, at["intended_target_name"]],
        rows[:, at["intended_target_chr"]],
        rows[:, at["intended_target_start"]],
        rows[:, at["intended_target_end"]],
    ):
        if name in coordinate_dict or name == "non-targeting":
            continue
        coordinate_dict[name] = [chrom, start, end]

    return coordinate_dict

def igv(mdata, gtf: str, method: Optional[Literal['sceptre', 'perturbo']] = None,
        results_key: str = 'test_results') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Generate bedpe and bedgraph data for generic or method-specific data"""
    # Set column names based on whether method is specified
    if not method:  # Generic columns
        log2_fc_col = "log2_fc"
        p_value_col = "p_value"
    else:  # Method-specific columns
        log2_fc_col = f"{method}_log2_fc"
        p_value_col = f"{method}_p_value"

    # Process coordinates
    coordinate_dict = process_coordinates(mdata)

    # Process GTF file
    df_gtf = read_gtf(gtf).to_pandas()
    gencode_df = df_gtf[['gene_id', 'gene_name']].copy()
    gencode_df['gene_id2'] = gencode_df['gene_id'].str.split('.').str[0]
    gencode_df = (
        gencode_df
        .rename(columns={'gene_name': 'gtf_gene_name'})
        .drop_duplicates(subset=['gene_id2'], keep='first')
    )

    # Initialize data structures
    bedpe = {}
    bedgraph = {}

    # Process test results. Only the four columns this function reads are
    # carried into the merge -- the table is twenty-odd columns of metrics wide,
    # and a 13.1M-row copy of the rest is pure overhead. A column that is absent
    # is left out so the dropna/merge below raise where they always raised.
    results = mdata.uns[results_key]
    wanted = [
        column
        for column in ('gene_id', 'intended_target_name', log2_fc_col, p_value_col)
        if column in results
    ]
    test_results = pd.DataFrame({k: results[k] for k in wanted})
    merged_df = test_results.merge(
        gencode_df[['gene_id2', 'gtf_gene_name']],
        left_on='gene_id',
        right_on='gene_id2',
        how='left'
    )

    # Filter out rows where required columns are missing
    merged_df = merged_df.dropna(subset=[log2_fc_col, p_value_col])

    target_names = merged_df["intended_target_name"].to_numpy(dtype=object)
    gene_ids = merged_df["gene_id"].to_numpy(dtype=object)
    gtf_gene_names = merged_df["gtf_gene_name"].to_numpy(dtype=object)

    # Object-dtype `==` is elementwise PyObject_RichCompare, so this is the same
    # verdict `row["intended_target_name"] == row["gtf_gene_name"]` gave.
    target_is_gene = target_names == gtf_gene_names
    target_placed = _contains(target_names, coordinate_dict)

    # PROMOTER interactions. The keys are only written when something matched,
    # so a run with no promoters still yields the column-less frame the
    # untouched defaultdict used to produce.
    promoter = target_is_gene & target_placed
    if promoter.any():
        coords = [coordinate_dict[name] for name in target_names[promoter]]
        bedgraph["chr"] = [entry[0] for entry in coords]
        bedgraph["start"] = [entry[1] for entry in coords]
        bedgraph["end"] = [entry[2] for entry in coords]
        bedgraph["p_value"] = _selected(merged_df[p_value_col], promoter)
        bedgraph["log2_fc"] = _selected(merged_df[log2_fc_col], promoter)

    # ENHANCER-GENE interactions
    enhancer = ~target_is_gene & target_placed & _contains(gene_ids, coordinate_dict)
    if enhancer.any():
        source = [coordinate_dict[name] for name in target_names[enhancer]]
        target = [coordinate_dict[gene_id] for gene_id in gene_ids[enhancer]]
        bedpe["chr1"] = [entry[0] for entry in source]
        bedpe["start1"] = [entry[1] for entry in source]
        bedpe["end1"] = [entry[2] for entry in source]
        bedpe["chr2"] = [entry[0] for entry in target]
        bedpe["start2"] = [entry[1] for entry in target]
        bedpe["end2"] = [entry[2] for entry in target]
        bedpe["p_value"] = _selected(merged_df[p_value_col], enhancer)
        bedpe["log2_fc"] = _selected(merged_df[log2_fc_col], enhancer)

    bedpe_df = pd.DataFrame(bedpe)
    bedgraph_df = pd.DataFrame(bedgraph)

    method_name = method if method else "Analysis"
    if bedpe_df.empty:
        print(f"Warning: {method_name} bedpe_df is empty.")
    if bedgraph_df.empty:
        print(f"Warning: {method_name} bedgraph_df is empty.")

    print(f"\n{method_name.capitalize()} statistics ({results_key}):")
    print(f"Number of enhancer-gene interactions: {len(bedpe_df)}")
    print(f"Number of promoter interactions: {len(bedgraph_df)}")

    return bedpe_df, bedgraph_df

def process_results_config(mdata, gtf: str, results_key: str, analysis_type: Optional[str] = None):
    """Process a single results configuration (either test_results or cis/trans_test_results)"""

    # Check if the key exists in mdata.uns
    if results_key not in mdata.uns:
        print(f"Warning: {results_key} not found in mdata.uns, skipping...")
        return

    print(f"\nProcessing {results_key}...")

    # Check available methods/columns
    results_df = pd.DataFrame(mdata.uns[results_key])
    cols = results_df.columns
    print(f"Available columns for {results_key}:", cols)

    # Create output directory
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)

    # Check for generic columns first
    if 'log2_fc' in cols and 'p_value' in cols:
        print(f"Using generic log2_fc and p_value columns for {results_key}")
        bedpe_df, bedgraph_df = igv(mdata, gtf, None, results_key)

        # Save generic files with analysis type prefix
        filename_prefix = f"{analysis_type}_" if analysis_type else ""
        bedpe_path = os.path.join(output_dir, f"{filename_prefix}evaluation.bedpe")
        bedgraph_path = os.path.join(output_dir, f"{filename_prefix}evaluation.bedgraph")

        bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
        bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)

        print(f"\nFiles saved for {results_key}:")
        print(f"bedpe file: {bedpe_path}")
        print(f"bedgraph file: {bedgraph_path}")
    else:
        # Check for method-specific columns
        available_methods = []
        if 'sceptre_log2_fc' in cols and not results_df['sceptre_log2_fc'].isna().all():
            available_methods.append('sceptre')
        if 'perturbo_log2_fc' in cols and not results_df['perturbo_log2_fc'].isna().all():
            available_methods.append('perturbo')

        print(f"Available methods for {results_key}: {available_methods}")

        if not available_methods:
            print(f"No methods with valid data found for {results_key}")
            return

        # Process data for available methods
        for method in available_methods:
            bedpe_df, bedgraph_df = igv(mdata, gtf, method, results_key)

            # Generate outputs with analysis type prefix
            filename_prefix = f"{analysis_type}_" if analysis_type else ""
            bedpe_path = os.path.join(output_dir, f"{filename_prefix}{method}.bedpe")
            bedgraph_path = os.path.join(output_dir, f"{filename_prefix}{method}.bedgraph")

            # Save files
            bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
            bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)

            print(f"\n{method.capitalize()} files saved for {results_key}:")
            print(f"bedpe file: {bedpe_path}")
            print(f"bedgraph file: {bedgraph_path}")

if __name__ == "__main__":
    print("Starting program...")
    parser = argparse.ArgumentParser(description="Process MuData and generate bedpe and bedgraph files")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--gtf", type=str, required=True, help="Path to the GTF file")
    parser.add_argument("--results_key", type=str, default="test_results",
                      help="Key for test results in mdata.uns")
    parser.add_argument("--default", action="store_true",
                      help="Process MuData with local- and global-analysis per-element results")

    args = parser.parse_args()

    print("Loading MuData file...")
    # backed="r" avoids loading gene/guide .X into memory; this script only
    # reads .var and .uns[results_key].
    mdata = mu.read(args.mdata_path, backed="r")

    # Nextflow declares this directory as the process output. Preserve an
    # explicit empty output when the selected inference mode has no plottable
    # result table instead of turning a successful no-op into a missing-output
    # task failure.
    os.makedirs("evaluation_output", exist_ok=True)

    # Determine which results to process based on --default flag
    if args.default:
        results_configs = [
            {"key": "local_analysis_per_element_results", "type": "local_analysis"},
            {"key": "global_analysis_per_element_results", "type": "global_analysis"}
        ]
    else:
        # Process single test_results
        results_configs = [
            {"key": args.results_key, "type": None}
        ]

    # Process each results configuration
    for config in results_configs:
        process_results_config(mdata, args.gtf, config["key"], config["type"])

    print("\nAll processing completed successfully")
