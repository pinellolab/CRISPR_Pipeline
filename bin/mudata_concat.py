#!/usr/bin/env python
"""Aggregate the guide-assigned MuData and decide, once, which cells are analysed.

This step runs immediately after guide assignment and is the only place the
pipeline chooses a cell population, so everything downstream -- both inference
methods, every covariate, both CRT pools -- inherits the same one. Before that
was true, SCEPTRE ran on ``bin/prepare_inference.py``'s cis subset (cells
carrying a tested or control guide) while PerTurbo fitted on every barcode in
the object, cells with no assigned guide included.

Dropping those cells would also destroy the guide-assignment rate, which is a
QC metric the dashboard reports: recomputed from the filtered object it reads
100% by construction, and a run whose guide assignment failed would look
perfect. So the counts as they stood before the filter are recorded in
``.uns`` under the names in ``ASSIGNED_GUIDE_FILTER_KEYS``, and the QC path
reads those rather than recounting.
"""
import argparse
import math
import os

import numpy as np
import pandas as pd
from scipy.sparse import issparse


ASSIGNMENT_LAYER = "guide_assignment"

# Recorded in the MuData's top-level .uns by
# ``filter_cells_without_assigned_guide``, whether or not the filter is applied,
# so the QC path reads one set of names either way. ``write_uns_patch`` copies
# .uns forward, and ``collapse_guides`` copies it too, so these survive into the
# published inference MuData.
ASSIGNED_GUIDE_FILTER_KEYS = (
    "assigned_guide_filter_applied",
    "assigned_guide_filter_source",
    "n_cells_before_assigned_guide_filter",
    "n_cells_after_assigned_guide_filter",
    "n_cells_with_assigned_guide",
    "n_cells_without_assigned_guide",
    "frac_cells_with_assigned_guide",
)


def parse_bool(value):
    """Accept Nextflow's rendered ``true``/``false`` as well as a real bool."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "t", "yes", "y", "1"}:
        return True
    if text in {"false", "f", "no", "n", "0"}:
        return False
    raise ValueError(f"Expected a boolean, got {value!r}.")


def assigned_guides_per_cell(mdata):
    """Assigned guides per cell, from the binarized guide-assignment layer.

    Returns the per-cell counts and the name of the matrix they came from. The
    layer is the guide *calls*; ``guide.X`` is raw guide UMIs and is only a
    fallback, because a cell with a stray guide read but no call would then
    count as assigned. Say so when that happens rather than silently switching.
    """
    guide = mdata.mod["guide"]
    if ASSIGNMENT_LAYER in guide.layers:
        matrix = guide.layers[ASSIGNMENT_LAYER]
        source = f"guide.layers['{ASSIGNMENT_LAYER}']"
    else:
        matrix = guide.X
        source = "guide.X"
        print(
            f"WARNING: the guide modality has no '{ASSIGNMENT_LAYER}' layer; counting "
            "assigned guides from raw guide UMIs in guide.X instead. A cell with guide "
            "reads but no call will be counted as assigned."
        )
    binarized = matrix > 0
    if issparse(binarized):
        counts = np.asarray(binarized.sum(axis=1)).ravel()
    else:
        counts = np.asarray(np.asarray(binarized).sum(axis=1)).ravel()
    return counts, source


def filter_cells_without_assigned_guide(mdata, require_assigned_guide=True):
    """Drop cells carrying no assigned guide, and record the pre-filter counts.

    A cell with exactly one assigned guide is kept; only ``guides_per_cell == 0``
    is dropped. The recorded counts describe the population as it stood here,
    which is the only point at which the guide-assignment rate is measurable,
    and they are written whether or not the filter is applied.
    """
    import mudata as md

    guides_per_cell, source = assigned_guides_per_cell(mdata)
    keep = guides_per_cell >= 1
    n_before = int(mdata.n_obs)
    n_with_guide = int(keep.sum())
    n_without_guide = n_before - n_with_guide
    n_after = n_with_guide if require_assigned_guide else n_before

    mdata.uns["assigned_guide_filter_applied"] = bool(require_assigned_guide)
    mdata.uns["assigned_guide_filter_source"] = source
    mdata.uns["n_cells_before_assigned_guide_filter"] = n_before
    mdata.uns["n_cells_after_assigned_guide_filter"] = n_after
    mdata.uns["n_cells_with_assigned_guide"] = n_with_guide
    mdata.uns["n_cells_without_assigned_guide"] = n_without_guide
    mdata.uns["frac_cells_with_assigned_guide"] = (
        float(n_with_guide) / n_before if n_before else 0.0
    )

    print(
        f"{n_with_guide} of {n_before} cells carry at least one assigned guide "
        f"({n_without_guide} carry none); counted from {source}"
    )

    if not require_assigned_guide:
        print(
            "QC_require_assigned_guide is off: keeping cells with no assigned guide. "
            "Both inference methods will fit on them."
        )
        return mdata

    if n_with_guide == 0:
        raise ValueError(
            "No cell carries an assigned guide, so the assigned-guide filter would "
            "empty the object. Check guide assignment, or set "
            "QC_require_assigned_guide = false to keep every cell."
        )

    if n_without_guide == 0:
        print("No cell to drop; every cell carries an assigned guide.")
        return mdata

    print(f"Dropping {n_without_guide} cells with no assigned guide")
    return mdata[keep].copy()


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

def apply_cell_and_gene_filters(mdata, min_cells_fraction, require_assigned_guide=True):
    """The cell filter, the per-cell depth totals, then the gene filter, in that order.

    Cells first: the gene threshold is a fraction of the *retained* cells, so a
    gene has to be supported among the cells that are actually analysed. On a
    screen with a large unassigned population the two orders disagree, and
    measuring gene support over cells nothing is tested in is the wrong one.

    The inference covariates are derived here because this is the first step
    that sees the final cell set: the depth totals are pinned before the gene
    filter (they must cover every gene), and the conditioned columns for both
    methods are written after it, so everything downstream inherits one set of
    values.
    """
    from inference_covariates import (
        ensure_gene_depth_totals,
        ensure_guide_umi_totals,
        materialize_shared_covariates,
    )

    mdata = filter_cells_without_assigned_guide(
        mdata, require_assigned_guide=require_assigned_guide
    )
    ensure_guide_umi_totals(mdata)
    ensure_gene_depth_totals(mdata)
    mdata = filter_genes_by_cells(mdata, min_cells_fraction)
    materialize_shared_covariates(mdata)
    return mdata


def concat_mudatas(
    input_files, output_file, min_cells_fraction=0.05, require_assigned_guide=True
):
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
        single_mdata = apply_cell_and_gene_filters(
            single_mdata, min_cells_fraction, require_assigned_guide
        )
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


    print('filtering cells and genes')
    combined_mdata = apply_cell_and_gene_filters(
        combined_mdata, min_cells_fraction, require_assigned_guide
    )


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
    parser.add_argument(
        "--require-assigned-guide",
        dest="require_assigned_guide",
        type=parse_bool,
        nargs="?",
        const=True,
        default=True,
        help=(
            "Keep only cells with at least one assigned guide (QC_require_assigned_guide). "
            "This is the pipeline's single cell filter, inherited by both inference "
            "methods. The pre-filter counts are recorded in .uns either way."
        ),
    )
    args = parser.parse_args()

    concat_mudatas(
        args.input,
        args.output,
        args.gene_filter,
        args.require_assigned_guide,
    )

if __name__ == "__main__":
    main()
