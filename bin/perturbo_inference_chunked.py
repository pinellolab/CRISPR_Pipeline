#!/usr/bin/env python
"""
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

from chunk_mudata import chunk_mudata
from perturbo_inference import resolve_efficiency_mode, resolve_num_workers


def get_gene_count(mdata_input_fp, gene_modality_name):
    mdata = md.read_h5mu(mdata_input_fp, backed="r")
    return mdata[gene_modality_name].n_vars


def should_chunk(n_genes, chunk_size):
    return chunk_size > 0 and n_genes > chunk_size


def build_perturbo_command(
    perturbo_script,
    mdata_input_fp,
    results_tsv_fp,
    mdata_output_fp=None,
    fit_guide_efficacy=True,
    efficiency_mode="scaled",
    accelerator="gpu",
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
    resolved_num_workers = resolve_num_workers(num_workers)
    resolved_efficiency_mode = resolve_efficiency_mode(efficiency_mode)
    cmd = [
        sys.executable,
        str(perturbo_script),
        mdata_input_fp,
        results_tsv_fp,
        "--fit_guide_efficacy",
        str(fit_guide_efficacy),
        "--efficiency_mode",
        resolved_efficiency_mode,
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
        "--inference_type",
        inference_type,
        "--num_workers",
        str(resolved_num_workers),
    ]
    if mdata_output_fp:
        cmd.extend(["--mdata_output_fp", mdata_output_fp])
    if test_all_pairs:
        cmd.append("--test_all_pairs")
    return cmd


def run_command(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, end="" if result.stderr.endswith("\n") else "\n")
        raise RuntimeError(f"Command failed with return code {result.returncode}")


def combine_chunk_results(result_files, output_path):
    dataframes = [pd.read_csv(result_file, sep="\t") for result_file in result_files]
    combined = pd.concat(dataframes, ignore_index=True) if dataframes else pd.DataFrame()

    if output_path.endswith(".gz"):
        combined.to_csv(output_path, index=False, sep="\t", compression="gzip")
    else:
        combined.to_csv(output_path, index=False, sep="\t")

    return combined


def write_results_mudata(base_mdata_path, output_path, inference_type, results_df):
    mdata = md.read(base_mdata_path)
    mdata.uns[f"per_{inference_type}_results"] = results_df
    mdata.write(output_path, compression="gzip")


def run_perturbo_chunked(
    mdata_input_fp,
    results_tsv_fp,
    mdata_output_fp=None,
    chunk_size=8000,
    fit_guide_efficacy=True,
    efficiency_mode="scaled",
    accelerator="gpu",
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
    script_dir = Path(__file__).resolve().parent
    perturbo_script = script_dir / "perturbo_inference.py"

    n_genes = get_gene_count(mdata_input_fp, gene_modality_name)
    if not should_chunk(n_genes, chunk_size):
        if chunk_size <= 0:
            print(
                f"Chunking disabled (chunk_size={chunk_size}); running standard PerTurbo inference."
            )
        else:
            print(
                f"Dataset has {n_genes} genes, which does not exceed chunk size {chunk_size}; running standard PerTurbo inference."
            )
        run_command(
            build_perturbo_command(
                perturbo_script=perturbo_script,
                mdata_input_fp=mdata_input_fp,
                results_tsv_fp=results_tsv_fp,
                mdata_output_fp=mdata_output_fp,
                fit_guide_efficacy=fit_guide_efficacy,
                efficiency_mode=efficiency_mode,
                accelerator=accelerator,
                batch_size=batch_size,
                early_stopping=early_stopping,
                early_stopping_patience=early_stopping_patience,
                lr=lr,
                num_epochs=num_epochs,
                gene_modality_name=gene_modality_name,
                guide_modality_name=guide_modality_name,
                inference_type=inference_type,
                test_all_pairs=test_all_pairs,
                num_workers=num_workers,
            )
        )
        return

    temp_dir = tempfile.mkdtemp(prefix="perturbo_chunks_")
    print(
        f"Chunking {n_genes} genes with max chunk size {chunk_size} into balanced subsets in {temp_dir}"
    )

    try:
        chunk_files = chunk_mudata(
            mudata_file=mdata_input_fp,
            output_dir=temp_dir,
            chunk_size=chunk_size,
            test_all_pairs=test_all_pairs,
            output_prefix="chunk",
        )

        result_files = []
        for chunk_file in chunk_files:
            chunk_result = os.path.splitext(chunk_file)[0] + ".tsv.gz"
            run_command(
                build_perturbo_command(
                    perturbo_script=perturbo_script,
                    mdata_input_fp=chunk_file,
                    results_tsv_fp=chunk_result,
                    fit_guide_efficacy=fit_guide_efficacy,
                    efficiency_mode=efficiency_mode,
                    accelerator=accelerator,
                    batch_size=batch_size,
                    early_stopping=early_stopping,
                    early_stopping_patience=early_stopping_patience,
                    lr=lr,
                    num_epochs=num_epochs,
                    gene_modality_name=gene_modality_name,
                    guide_modality_name=guide_modality_name,
                    inference_type=inference_type,
                    test_all_pairs=test_all_pairs,
                    num_workers=num_workers,
                )
            )
            result_files.append(chunk_result)

        if not result_files:
            raise RuntimeError("Chunked PerTurbo did not produce any result files")

        combined_results = combine_chunk_results(result_files, results_tsv_fp)
        if mdata_output_fp:
            write_results_mudata(
                base_mdata_path=mdata_input_fp,
                output_path=mdata_output_fp,
                inference_type=inference_type,
                results_df=combined_results,
            )
        print(f"Wrote combined chunked results to {results_tsv_fp}")
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
        choices=["scaled"],
        default="scaled",
        help="Efficiency mode for the model. Only 'scaled' is supported.",
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
