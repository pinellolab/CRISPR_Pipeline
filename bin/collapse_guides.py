#!/usr/bin/env python
import argparse
import mudata as md


def collapse_guides(
    mdata_input_fp,
    mdata_output_fp,
    # Add additional parameters here as needed
):
    """
    Collapse guides in MuData object.

    Parameters:
    -----------
    mdata_input_fp : str
        Path to input MuData file
    mdata_output_fp : str
        Path to output MuData file
    """
    mdata = md.read(mdata_input_fp)

    # TODO: Implement guide collapsing logic here

    mdata.write(mdata_output_fp)
    return mdata


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
