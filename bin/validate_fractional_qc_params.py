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

    for key, default in (
        ("QC_min_counts_per_cell", 500),
        ("QC_MAD_total_counts", 5),
        ("QC_MAD_n_genes", 5),
    ):
        try:
            numeric = float(params.get(key, default))
        except (TypeError, ValueError) as error:
            raise ValueError(f"{key} must be numeric.") from error
        if numeric < 0:
            raise ValueError(f"{key} must be non-negative.")
    pct_mito = float(params.get("QC_pct_mito", 15))
    if not 0 <= pct_mito <= 100:
        raise ValueError("QC_pct_mito must be in [0, 100].")
    scrublet_profile = params.get("SCRUBLET_assay_type", "droplet")
    if scrublet_profile not in {"droplet", "cc-perturb-seq"}:
        raise ValueError("SCRUBLET_assay_type must be droplet or cc-perturb-seq.")
    scrublet_rate = params.get("SCRUBLET_expected_doublet_rate")
    if scrublet_rate is not None and not 0 < float(scrublet_rate) < 1:
        raise ValueError("SCRUBLET_expected_doublet_rate must be null or in (0, 1).")
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
