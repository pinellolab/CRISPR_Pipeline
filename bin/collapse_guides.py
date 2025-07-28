#!/usr/bin/env python
import argparse
from re import M
import mudata as md
import pandas as pd
import numpy as np
from scipy.sparse import issparse


def collapse_guides(
    mdata_input_fp,
    mdata_output_fp,
    max_elements_per_cell=1,
    max_guides_per_cell=2,
    min_guides_per_cell=2,
):
    """
    Collapse guides in MuData object.

    This function identifies cells with exactly 2 guides targeting the same element
    and collapses them into a single guide entry with combined guide_ids.

    Parameters:
    -----------
    mdata_input_fp : str
        Path to input MuData file
    mdata_output_fp : str
        Path to output MuData file
    """
    mdata = md.read(mdata_input_fp)
    guide_adata = mdata["guide"]

    # Get guides per cell from guide assignment layer
    guides_per_cell = mdata["guide"].layers["guide_assignment"].sum(axis=1)
    if issparse(guides_per_cell):
        guides_per_cell = guides_per_cell.A1
    elif hasattr(guides_per_cell, "values"):
        guides_per_cell = guides_per_cell.values.flatten()
    else:
        guides_per_cell = np.asarray(guides_per_cell).flatten()

    # Create a guide by element matrix
    guide_by_element = pd.get_dummies(guide_adata.var["intended_target_name"])
    cell_by_element = (mdata["guide"].layers["guide_assignment"] @ guide_by_element) > 0

    # Get elements per cell as 1D array
    elements_per_cell = cell_by_element.sum(axis=1)
    if issparse(elements_per_cell):
        elements_per_cell = elements_per_cell.A1
    elif hasattr(elements_per_cell, "values"):
        elements_per_cell = elements_per_cell.values.flatten()
    else:
        elements_per_cell = np.asarray(elements_per_cell).flatten()

    # Filter for cells with no more than max_elements_per_cell elements and specified number of guides
    element_filter = elements_per_cell <= max_elements_per_cell
    guide_filter = (guides_per_cell >= min_guides_per_cell) & (
        guides_per_cell <= max_guides_per_cell
    )
    filtered_cells = element_filter & guide_filter
    print(
        f"Found {filtered_cells.sum()} cells with {min_guides_per_cell} <= guides <= {max_guides_per_cell} and elements <= {max_elements_per_cell}."
    )

    guide_adata_subset = guide_adata[filtered_cells, :]

    # Get nonzero indices for assigned guides
    if issparse(guide_adata_subset.layers["guide_assignment"]):
        cell_idx, guide_idx = guide_adata_subset.layers["guide_assignment"].nonzero()
    else:
        cell_idx, guide_idx = np.where(guide_adata_subset.layers["guide_assignment"])

    # Create dataframe with cell barcodes and guide information
    bc_df = pd.DataFrame({"cell_bc": guide_adata_subset.obs_names[cell_idx]})
    var_df = guide_adata_subset.var.iloc[guide_idx].reset_index(drop=True)
    bc_var_df = pd.concat([bc_df, var_df], axis=1)

    # Group by cell and target information, then collapse guide_ids
    grouping_cols = [
        "cell_bc",
        "intended_target_name",
    ]

    bc_var_df_grouped = bc_var_df.sort_values("guide_id").groupby(
        grouping_cols, observed=True
    )

    guide_spacer_grouped = bc_var_df_grouped.agg(
        {
            "guide_id": lambda x: "|".join(x.tolist()),
            "spacer": lambda x: "|".join(x.tolist()),
        }
    )

    # Keep first value for other columns
    rest_cols = [
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
        "targeting",
        "type",
        "pam",
        "gene_name",
    ]
    rest_grouped = bc_var_df_grouped[rest_cols].agg("first")

    # Combine collapsed data
    cell_guide_deduped = (
        pd.concat([guide_spacer_grouped, rest_grouped], axis=1)
        .reset_index()
        .set_index("cell_bc")
        .assign(
            guide_chr=lambda x: x["intended_target_chr"],
            guide_start=lambda x: x["intended_target_start"],
            guide_end=lambda x: x["intended_target_end"],
        )
    )

    # Create new guide matrix with collapsed guides
    guide_dummies = pd.get_dummies(cell_guide_deduped["guide_id"], sparse=True)
    print(
        f"Collapsed {len(bc_var_df)} guide assignments into {len(guide_dummies.columns)} unique guide combinations"
    )

    var_new = cell_guide_deduped.groupby("guide_id").first()

    # Create new AnnData object for collapsed guides
    collapsed_guide_adata = md.AnnData(
        guide_dummies,
        obs=mdata["guide"].obs.loc[guide_dummies.index],
        var=var_new.loc[guide_dummies.columns],
    )

    # Add guide metadata to the new var dataframe

    # for guide_id in collapsed_guide_adata.var_names:
    #     # Find the first occurrence of this guide_id (or combination) in the deduped data
    #     mask = cell_guide_deduped["guide_id"] == guide_id
    #     if mask.any():
    #         guide_info = cell_guide_deduped[mask].iloc[0]
    #         for col in rest_cols + grouping_cols[1:]:  # Skip cell_bc
    #             if col in guide_info:
    #                 val = guide_info[col]
    #                 # Convert to string if it's a categorical or non-string object
    #                 if hasattr(val, "name"):  # pandas categorical
    #                     val = str(val)
    #                 elif not isinstance(val, (str, int, float, bool, np.number)):
    #                     val = str(val)
    #                 collapsed_guide_adata.var.loc[guide_id, col] = val

    # Create guide assignment layer for collapsed guides

    collapsed_guide_adata.layers["guide_assignment"] = collapsed_guide_adata.X.copy()

    # include the raw guide counts with the collapsed guide adata
    collapsed_guide_adata.raw = mdata["guide"][
        collapsed_guide_adata.obs_names, :
    ].copy()

    # Create new MuData object with collapsed guides
    mdata_dict = {}
    for modality in mdata.mod.keys():
        if modality == "guide":
            mdata_dict[modality] = collapsed_guide_adata
        else:
            mdata_dict[modality] = mdata[modality][
                collapsed_guide_adata.obs_names, :
            ].copy()

    mdata_collapsed = md.MuData(mdata_dict)

    mdata_collapsed.write(mdata_output_fp)
    return mdata_collapsed


def main():
    parser = argparse.ArgumentParser(description="Collapse guides in MuData object")
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument("mdata_output_fp", type=str, help="Output file path for MuData")

    # Add additional arguments as needed
    # parser.add_argument(
    #     "--example_param",
    #     type=str,
    #     default="default_value",
    #     help="Example parameter description"
    # )

    args = parser.parse_args()
    collapse_guides(
        args.mdata_input_fp,
        args.mdata_output_fp,
    )


if __name__ == "__main__":
    main()
