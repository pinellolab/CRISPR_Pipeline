#!/usr/bin/env python
import argparse
import mudata as md
import numpy as np


def create_nontargeting_intended_targets(
    mdata_input_fp,
    mdata_output_fp,
    strategy="same",  # strategy for creating NT elements: "median" or "same"
):
    """
    Fill intended target names for nontargeting/safe-targeting/negative control guides in MuData object.

    This function assigns intended target names to control guides (non-targeting,
    safe-targeting, negative control) based on the specified strategy.

    Parameters:
    -----------
    mdata_input_fp : str
        Path to input MuData file
    mdata_output_fp : str
        Path to output MuData file
    strategy : str
        Strategy for assigning targets ("median" or "same")
        If median, creates new NT elements from the NT guides with the same median guides per element as the targeting guides.
        If same, assigns all NT guides to a single "non-targeting" element.
    """
    mdata = md.read(mdata_input_fp)
    gdata = mdata["guide"]

    # Identify control guides
    control_guides = (~gdata.var["targeting"]) | (
        gdata.var["type"].isin(["safe-targeting", "non-targeting", "negative control"])
    )

    control_guide_var = gdata.var[control_guides]
    targeting_guide_var = gdata.var[~control_guides]

    print(
        f"Found {control_guide_var.shape[0]} control guides and {targeting_guide_var.shape[0]} targeting guides"
    )

    # Calculate guides per element for targeting guides
    guides_per_element = (
        targeting_guide_var.groupby(["intended_target_name"], observed=True)
        .size()
        .sort_values()
    )

    print("Guides per element distribution:")
    print(guides_per_element.value_counts().sort_index())

    # Generate intended targets based on strategy
    if strategy == "median":
        median_guides_per_element = int(guides_per_element.median())
        total_nt_guides = control_guide_var.shape[0]
        total_nt_elements = int(np.ceil(total_nt_guides / median_guides_per_element))

        print(f"Using median strategy: {median_guides_per_element} guides per element")
        print(f"Creating {total_nt_elements} non-targeting elements")

        intended_targets = [
            f"non-targeting|{i + 1}"
            for i in range(total_nt_elements)
            for _ in range(median_guides_per_element)
        ][:total_nt_guides]  # Trim to exact number needed

    elif strategy == "same":
        print("Using same strategy: all control guides assigned to 'non-targeting'")
        intended_targets = ["non-targeting"] * control_guide_var.shape[0]
    else:
        raise ValueError(f"Unknown strategy: {strategy}. Must be 'median' or 'same'")

    # Update intended target names
    gdata.var["intended_target_name"] = gdata.var["intended_target_name"].astype("str")
    gdata.var.loc[control_guides, "intended_target_name"] = intended_targets
    gdata.var["intended_target_name"] = gdata.var["intended_target_name"].astype(
        "category"
    )

    # Update the MuData object
    mdata.update()

    print(f"Updated {len(intended_targets)} control guides with intended target names")

    mdata.write(mdata_output_fp)
    return mdata


def main():
    parser = argparse.ArgumentParser(
        description="Fill missing intended target names for control guides"
    )
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument("mdata_output_fp", type=str, help="Output file path for MuData")
    parser.add_argument(
        "--strategy",
        type=str,
        choices=["median", "same"],
        default="median",
        help="Strategy for assigning targets to control guides (default: median)",
    )

    args = parser.parse_args()
    create_nontargeting_intended_targets(
        args.mdata_input_fp,
        args.mdata_output_fp,
        strategy=args.strategy,
    )


if __name__ == "__main__":
    main()
