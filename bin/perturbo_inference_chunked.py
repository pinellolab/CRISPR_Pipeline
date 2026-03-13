#!/usr/bin/env python
"""
Run PerTurbo inference on balanced gene chunks and concatenate TSV outputs.
Run PerTurbo inference on balanced gene chunks and concatenate TSV outputs.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import mudata as md
import pandas as pd
import mudata as md


def run_perturbo_chunked(
    mdata_input_fp,
    results_tsv_fp,
    mdata_output_fp=None,
    fit_guide_efficacy=True,
    efficiency_mode="undecided",
    accelerator="auto",
    batch_size=4096,
    early_stopping=False,
    early_stopping_patience=5,
    lr=0.01,
    num_epochs=100,
    gene_modality_name="gene",
    guide_modality_name="guide",
    inference_type="element",
    test_all_pairs=False,
    num_workers=None,
):
    """
    Run PerTurbo inference on chunks of genes by calling external scripts.

    This function:
    1. Uses chunk_mudata.py to create gene chunks
    2. Calls perturbo_inference.py on each chunk via subprocess
    3. Combines the results and saves them back to the original MuData

    Parameters:
    -----------
    mdata_input_fp : str
        Path to input MuData file
    mdata_output_fp : str
        Path to output MuData file
    chunk_size : int
        Number of genes per chunk (default: 8000)
    Other parameters match those in perturbo_inference.py
    """

    print(f"Starting chunked PerTurbo inference on {mdata_input_fp}...")
    mdata = md.read_h5mu(mdata_input_fp, backed="r")
    if mdata[gene_modality_name].n_vars <= chunk_size:
        print(
            "Number of genes is less than or equal to chunk size; running standard PerTurbo inference instead."
        )
        perturbo_script = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "perturbo_inference.py"
        )
        perturbo_cmd = [
            "python",
            perturbo_script,
            mdata_input_fp,
            mdata_output_fp,
            "--fit_guide_efficacy",
            str(fit_guide_efficacy),
            "--efficiency_mode",
            efficiency_mode,
            "--accelerator",
            accelerator,
            "--batch_size",
            str(batch_size),
            "--early_stopping",
            str(early_stopping),
            "--early_stopping_patience",
            str(early_stopping_patience),
            "--lr",
            str(lr),
            "--num_epochs",
            str(num_epochs),
            "--gene_modality_name",
            gene_modality_name,
            "--guide_modality_name",
            guide_modality_name,
            "--test_control_guides",
            str(test_control_guides),
            "--num_workers",
            str(num_workers),
        ]

        if test_all_pairs:
            perturbo_cmd.append("--test_all_pairs")

        print(f"Running: {' '.join(perturbo_cmd)}")
        subprocess.run(perturbo_cmd)
        return

    temp_dir = tempfile.mkdtemp(prefix="perturbo_chunks_")
    print(f"Using temporary directory: {temp_dir}")

    try:
        # Step 1: Create chunks using chunk_mudata.py
        print(f"\nStep 1: Chunking MuData into {chunk_size} genes per chunk...")
        chunk_cmd = [
            "python",
            chunk_script,
            mdata_input_fp,
            temp_dir,
            "--chunk-size",
            str(chunk_size),
            "--output-prefix",
            "chunk",
        ]

        if test_all_pairs:
            chunk_cmd.append(
                "--test-all-pairs"
            )  # Add flag when we want to test all pairs

        print(f"Running: {' '.join(chunk_cmd)}")
        result = subprocess.run(chunk_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error in chunking: {result.stderr}")
            raise RuntimeError(f"Chunking failed with return code {result.returncode}")

        print("Chunking completed successfully")

        # Step 2: Find all chunk files
        chunk_files = []
        for filename in os.listdir(temp_dir):
            if filename.startswith("chunk.") and filename.endswith(".h5mu"):
                chunk_files.append(os.path.join(temp_dir, filename))

        chunk_files.sort()  # Ensure consistent order
        print(f"Found {len(chunk_files)} chunk files")

        if not chunk_files:
            raise RuntimeError("No chunk files were created")

        # Step 3: Run PerTurbo on each chunk
        processed_chunks = []

        for i, chunk_file in enumerate(chunk_files):
            print(
                f"\nStep 3.{i + 1}: Processing chunk {i + 1}/{len(chunk_files)}: {os.path.basename(chunk_file)}"
            )

            # Create output file for this chunk
            chunk_output = chunk_file.replace(".h5mu", "_results.h5mu")

            # Build command for perturbo_inference.py
            perturbo_cmd = [
                "python",
                perturbo_script,
                chunk_file,
                chunk_output,
                "--fit_guide_efficacy",
                str(fit_guide_efficacy),
                "--efficiency_mode",
                efficiency_mode,
                "--accelerator",
                accelerator,
                "--batch_size",
                str(batch_size),
                "--early_stopping",
                str(early_stopping),
                "--early_stopping_patience",
                str(early_stopping_patience),
                "--lr",
                str(lr),
                "--num_epochs",
                str(num_epochs),
                "--gene_modality_name",
                gene_modality_name,
                "--guide_modality_name",
                guide_modality_name,
                "--test_control_guides",
                str(test_control_guides),
                "--num_workers",
                str(num_workers),
            ]

            if test_all_pairs:
                perturbo_cmd.append("--test_all_pairs")

            print(f"Running: {' '.join(perturbo_cmd)}")
            result = subprocess.run(perturbo_cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"Error in chunk {i + 1}: {result.stderr}")
                print(f"Skipping chunk {i + 1}")
                continue

            processed_chunks.append(chunk_output)
            print(f"Successfully processed chunk {i + 1}")

        if not processed_chunks:
            raise RuntimeError("No chunks were successfully processed")

        # Step 4: Combine results
        print(
            f"\nStep 4: Combining results from {len(processed_chunks)} processed chunks..."
        )
        combine_chunk_results(processed_chunks, mdata_input_fp, mdata_output_fp)

        print(f"Results saved to {mdata_output_fp}")

    finally:
        shutil.rmtree(temp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Run chunked PerTurbo analysis on MuData",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument(
        "results_tsv_fp",
        type=str,
        help="Output TSV file path for concatenated PerTurbo results",
    )
    parser.add_argument(
        "--mdata_output_fp",
        type=str,
        default=None,
        help="Optional output file path for MuData; if omitted the MuData will not be written",
    )
    parser.add_argument(
        "--chunk_size",
        "-c",
        type=int,
        default=8000,
        help="Maximum genes per chunk; values <= 0 disable chunking",
    )
    parser.add_argument(
        "--fit_guide_efficacy",
        type=bool,
        default=True,
        help="Whether to fit guide efficacy",
    )
    parser.add_argument(
        "--efficiency_mode",
        type=str,
        choices=["undecided", "low", "high"],
        default="undecided",
        help="Efficiency mode for the model: 'undecided'/'auto' (auto-detect), 'low' (mixture), 'high' (scaled)",
    )
    parser.add_argument(
        "--accelerator",
        type=str,
        choices=["auto", "gpu", "cpu"],
        default="gpu",
        help="Accelerator to use for training",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=4096,
        help="Batch size for training",
    )
    parser.add_argument(
        "--early_stopping",
        type=bool,
        default=False,
        help="Whether to use early stopping during training",
    )
    parser.add_argument(
        "--early_stopping_patience",
        type=int,
        default=5,
        help="Patience for early stopping",
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=0.01,
        help="Learning rate for training",
    )
    parser.add_argument(
        "--num_epochs",
        type=int,
        default=100,
        help="Maximum number of epochs for training",
    )
    parser.add_argument(
        "--gene_modality_name",
        type=str,
        default="gene",
        help="Name of the gene modality in the MuData object",
    )
    parser.add_argument(
        "--guide_modality_name",
        type=str,
        default="guide",
        help="Name of the guide modality in the MuData object",
    )
    parser.add_argument(
        "--inference_type",
        type=str,
        default="element",
        help="Unit to test for effects on each gene: 'guide' or 'element'",
    )
    parser.add_argument(
        "--inference_type",
        type=str,
        default="element",
        help="Unit to test for effects on each gene: 'guide' or 'element'",
    )
    parser.add_argument(
        "--test_all_pairs",
        action="store_true",
        help="Whether to test all pairs or only those in pairs_to_test",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=None,
        help="Number of workers for data loading (default: assigned CPUs minus one)",
    )
    args = parser.parse_args()

    run_perturbo_chunked(
        mdata_input_fp=args.mdata_input_fp,
        results_tsv_fp=args.results_tsv_fp,
        mdata_output_fp=args.mdata_output_fp,
        chunk_size=args.chunk_size,
        fit_guide_efficacy=args.fit_guide_efficacy,
        efficiency_mode=args.efficiency_mode,
        accelerator=args.accelerator,
        batch_size=args.batch_size,
        early_stopping=args.early_stopping,
        early_stopping_patience=args.early_stopping_patience,
        lr=args.lr,
        num_epochs=args.num_epochs,
        gene_modality_name=args.gene_modality_name,
        guide_modality_name=args.guide_modality_name,
        inference_type=args.inference_type,
        test_all_pairs=args.test_all_pairs,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
