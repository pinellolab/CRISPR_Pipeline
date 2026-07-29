#!/usr/bin/env python3
"""Fail fast when a run still uses the retired absolute gene-support setting."""

import argparse
import json


def validate_fractional_qc_params(params):
    if "QC_min_cells_per_gene" not in params:
        raise ValueError("Missing required parameter QC_min_cells_per_gene.")
    try:
        value = float(params["QC_min_cells_per_gene"])
    except (TypeError, ValueError) as error:
        raise ValueError("QC_min_cells_per_gene must be numeric.") from error
    if not 0 <= value < 1:
        raise ValueError(
            "QC_min_cells_per_gene must be a fraction in [0, 1); "
            f"received {params['QC_min_cells_per_gene']!r}."
        )

    tapseq_mode = params.get("TAPSEQ_QC_MODE", False)
    if not isinstance(tapseq_mode, bool):
        raise ValueError("TAPSEQ_QC_MODE must be true or false.")
    return value, tapseq_mode


def main():
    parser = argparse.ArgumentParser(
        description="Validate fractional gene-support parameters before Nextflow starts."
    )
    parser.add_argument("params_json")
    args = parser.parse_args()
    with open(args.params_json) as handle:
        params = json.load(handle)
    value, tapseq_mode = validate_fractional_qc_params(params)
    print(
        "Fractional QC parameters validated: "
        f"QC_min_cells_per_gene={value:g}, TAPSEQ_QC_MODE={str(tapseq_mode).lower()}"
    )


if __name__ == "__main__":
    main()
