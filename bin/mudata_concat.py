#!/usr/bin/env python
import argparse
import math
import os

import pandas as pd


def resolve_min_cells(n_obs, min_cells_fraction):
    """
    Resolve a fractional gene-support threshold.

    The threshold is fraction-only and retains the historical
    strict-greater-than behavior. A zero fraction keeps every gene detected
    in at least one cell.
    """
    min_cells_fraction = float(min_cells_fraction)
    if not 0 <= min_cells_fraction < 1:
        raise ValueError("Gene cell-support threshold must be a fraction in [0, 1).")
    return max(1, math.floor(n_obs * min_cells_fraction) + 1)


def filter_genes_by_cells(mdata, min_cells_fraction):
    """
    Filter genes by the minimum fraction of cells expressing them.
    """
    required_cells = resolve_min_cells(mdata['gene'].n_obs, min_cells_fraction)
    detected_cells = (mdata['gene'].X > 0).sum(0).A1
    index_filter = detected_cells >= required_cells
    print(
        f"Keeping {int(index_filter.sum())} of {len(index_filter)} genes "
        f"detected in at least {required_cells} cells"
    )
    mdata.mod['gene'] = mdata.mod['gene'][:, index_filter]
    return mdata


def preserve_source_guide_metadata(combined_guide_var, source_guide_var):
    """Restore source-only guide annotations after MuData concatenation.

    ``mudata.concat`` only guarantees the shared annotation schema. Library
    metadata such as an explicit ``element_id`` can therefore disappear even
    when every input has it. Those identifiers encode real paired constructs
    in dual-guide assays, so restore source-only columns without coercing the
    dtypes of annotations already handled by MuData.
    """
    combined = combined_guide_var.copy()
    source = source_guide_var.reindex(combined.index)

    for column in source.columns:
        if column not in combined.columns:
            combined[column] = source[column]

    return combined

def concat_mudatas(input_files, output_file, min_cells_fraction=0.05):
    """
    Concatenate multiple MuData files. If only one file is provided, it's copied to the output.
    """
    import mudata as md

    files = sorted(input_files, key=lambda x: os.path.basename(x))
    if not files:
        print(f"No files found: {input_files}")
        return

    print(f"Found {len(files)} files to concatenate")

    # Handle single file case
    if len(files) == 1:
        print(f"Only one file found. Copying {files[0]} to {output_file}")
        single_mdata = md.read(files[0])
        single_mdata = filter_genes_by_cells(single_mdata, min_cells_fraction)  # Filter genes based on minimum cells
        print(f"Saving MuData with {single_mdata.n_obs} cells to {output_file}")
        single_mdata.write(output_file)
        return

    # Handle multiple files case
    print("Concatenating all MuData objects...")
    mudatas = [md.read(ff) for ff in files]
    combined_mdata = md.concat(mudatas, merge='first', uns_merge='first', join='outer')

    # Keep assay/library annotations such as element_id. Collapse later uses
    # them to keep the two guides of an explicit control construct together.
    combined_mdata.mod['guide'].var = preserve_source_guide_metadata(
        combined_mdata.mod['guide'].var,
        mudatas[0].mod['guide'].var,
    )


    print ('filtering genes')
    combined_mdata = filter_genes_by_cells(combined_mdata, min_cells_fraction)  # Filter genes based on minimum cells


    print(f"Saving combined MuData with {combined_mdata.n_obs} cells to {output_file}")
    combined_mdata.write(output_file)

    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Concatenate MuData files")
    parser.add_argument("-i", "--input", dest="input", nargs="+", required=True, help="Input mudata files")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output file path")
    parser.add_argument(
        "-g",
        "--gene_filter",
        dest="gene_filter",
        type=float,
        default=0.05,
        help=(
            "Fraction of retained cells required to keep a gene. Must be in "
            "[0, 1); zero keeps every gene detected in at least one cell."
        ),
    )
    args = parser.parse_args()

    concat_mudatas(args.input, args.output, args.gene_filter)

if __name__ == "__main__":
    main()
