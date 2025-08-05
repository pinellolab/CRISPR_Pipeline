#!/usr/bin/env python
"""
Chunked PerTurbo inference script.

This script performs PerTurbo inference on chunks of genes to handle large datasets.
It adds covariates to the full dataset, chunks the data, runs inference on each chunk
independently, then combines the results and saves them back to the original MuData.
"""

import argparse
import os
from math import ceil
import pandas as pd
import numpy as np
import mudata as md
import perturbo
import scvi


def run_perturbo_chunked(
    mdata_input_fp,
    mdata_output_fp,
    chunk_size=4000,
    fit_guide_efficacy=True,
    efficiency_mode="undecided",
    accelerator="auto",
    batch_size=512,
    early_stopping=True,
    early_stopping_patience=5,
    lr=0.01,
    num_epochs=100,
    gene_modality_name="gene",
    guide_modality_name="guide",
    test_all_pairs=True,
    test_control_guides=True,
    num_workers=0,
):
    """
    Run PerTurbo inference on chunks of genes and combine results.

    Parameters:
    -----------
    mdata_input_fp : str
        Path to input MuData file
    mdata_output_fp : str
        Path to output MuData file
    chunk_size : int
        Number of genes per chunk (default: 4000)
    Other parameters match those in perturbo_inference.py
    """

    scvi.settings.seed = 0
    if num_workers > 0:
        scvi.settings.dl_num_workers = num_workers

    print(f"Loading MuData from {mdata_input_fp}...")
    mdata = md.read(mdata_input_fp)

    # Add covariates to the full dataset (as done in original perturbo_inference.py)
    print("Adding covariates to full dataset...")

    # Identify control guides
    control_guide_filter = (~mdata[guide_modality_name].var["targeting"]) | (
        mdata[guide_modality_name]
        .var["type"]
        .isin(["safe-targeting", "non-targeting", "negative control"])
    )
    if np.any(control_guide_filter):
        control_guides = (
            mdata[guide_modality_name].var_names[control_guide_filter].tolist()
        )
    else:
        control_guides = None

    # Add log-transformed covariates
    mdata[gene_modality_name].obs["log1p_total_guide_umis"] = np.log1p(
        mdata[guide_modality_name].obs["total_guide_umis"]
    )
    mdata[gene_modality_name].obs["log1p_total_guide_umis_centered"] = (
        mdata[gene_modality_name].obs["log1p_total_guide_umis"]
        - mdata[gene_modality_name].obs["log1p_total_guide_umis"].mean()
    )

    # Determine efficiency mode
    efficiency_mode = {"undecided": "auto", "low": "mixture", "high": "scaled"}[
        efficiency_mode
    ]

    if efficiency_mode == "auto":
        max_guides_per_cell = mdata[guide_modality_name].X.sum(axis=1).max()
        if max_guides_per_cell > 1:
            efficiency_mode = "scaled"
            print(
                "Using 'scaled' efficiency mode due to high MOI (max guides per cell > 1)."
            )
        else:
            efficiency_mode = "mixture"
            print(
                "Using 'mixture' efficiency mode due to low MOI (max guides per cell <= 1)."
            )

    # Create intended_targets matrix
    intended_targets_df = pd.get_dummies(
        mdata[guide_modality_name].var["intended_target_name"]
    ).astype(float)

    if not test_control_guides:
        if control_guides is not None:
            control_elements = intended_targets_df.loc[control_guides].sum(axis=0) > 0
            intended_targets_df = intended_targets_df.drop(
                columns=intended_targets_df.columns[control_elements]
            )

    mdata[guide_modality_name].varm["intended_targets"] = intended_targets_df
    mdata.uns["intended_target_names"] = intended_targets_df.columns.tolist()

    # Load pairs_to_test if needed
    pairs_to_test_df = None
    if not test_all_pairs:
        if "pairs_to_test" not in mdata.uns:
            raise ValueError(
                "pairs_to_test not found in mdata.uns when test_all_pairs=False"
            )

        pairs_to_test = mdata.uns["pairs_to_test"]
        if isinstance(pairs_to_test, dict):
            pairs_to_test_df = pd.DataFrame(pairs_to_test)
        else:
            pairs_to_test_df = pairs_to_test

    # Set up chunking
    n_genes = mdata[gene_modality_name].shape[1]
    n_chunks = ceil(n_genes / chunk_size)
    print(f"Splitting {n_genes} genes into {n_chunks} chunks of size {chunk_size}")

    all_results = []

    # Process each chunk
    for i in range(n_chunks):
        print(f"\n--- Processing chunk {i + 1}/{n_chunks} ---")

        # Define gene range for this chunk
        start = i * chunk_size
        end = min(start + chunk_size, n_genes)
        print(f"Genes {start} to {end - 1}")

        # Extract gene subset for this chunk
        gene_adata_subset = mdata[gene_modality_name][:, start:end]

        if not test_all_pairs:
            # Filter to genes and guides in this chunk
            chunk_genes = gene_adata_subset.var_names.tolist()
            chunk_pairs = pairs_to_test_df[
                pairs_to_test_df["gene_id"].isin(chunk_genes)
            ]

            if len(chunk_pairs) == 0:
                print(f"  No pairs to test in chunk {i + 1}, skipping...")
                continue

            chunk_guides = chunk_pairs["guide_id"].unique()
            print(
                f"  Found {len(chunk_guides)} unique guides for {len(chunk_genes)} genes"
            )

            # Filter guide data
            guide_mask = mdata[guide_modality_name].var["guide_id"].isin(chunk_guides)
            guide_adata_subset = mdata[guide_modality_name][:, guide_mask]

            # Create gene by element matrix for this chunk
            aggregated_df = (
                chunk_pairs[["gene_id", "intended_target_name"]]
                .drop_duplicates()
                .assign(value=1)
            )

            gene_by_element = (
                aggregated_df.pivot(
                    index="gene_id", columns="intended_target_name", values="value"
                )
                .reindex(
                    index=gene_adata_subset.var_names,
                    columns=mdata.uns["intended_target_names"],
                )
                .fillna(0)
            )

            # Subset to targeted cells only
            targeted_cells = guide_adata_subset.X.sum(axis=1) > 0
            if not targeted_cells.any():
                print(f"  No targeted cells found in chunk {i + 1}, skipping...")
                continue

            chunk_mdata = md.MuData(
                {
                    gene_modality_name: gene_adata_subset[targeted_cells],
                    guide_modality_name: guide_adata_subset[targeted_cells],
                }
            )

            # Add gene by element matrix to chunk
            chunk_mdata[gene_modality_name].varm["intended_targets"] = gene_by_element
            chunk_mdata.uns["pairs_to_test"] = chunk_pairs

        else:
            # Use all guide data for this gene chunk
            chunk_mdata = md.MuData(
                {
                    gene_modality_name: gene_adata_subset,
                    guide_modality_name: mdata[guide_modality_name].copy(),
                }
            )

        # Copy over necessary metadata
        chunk_mdata.uns["intended_target_names"] = mdata.uns["intended_target_names"]

        print(
            f"  Chunk shape: {gene_modality_name} {chunk_mdata[gene_modality_name].shape}, "
            f"{guide_modality_name} {chunk_mdata[guide_modality_name].shape}"
        )

        try:
            # Run PerTurbo on this chunk
            chunk_results = run_perturbo_on_chunk(
                chunk_mdata,
                efficiency_mode=efficiency_mode,
                fit_guide_efficacy=fit_guide_efficacy,
                accelerator=accelerator,
                batch_size=batch_size,
                early_stopping=early_stopping,
                early_stopping_patience=early_stopping_patience,
                lr=lr,
                num_epochs=num_epochs,
                gene_modality_name=gene_modality_name,
                guide_modality_name=guide_modality_name,
                test_all_pairs=test_all_pairs,
                control_guides=control_guides,
            )

            if chunk_results is not None:
                all_results.append(chunk_results)
                print(
                    f"  Successfully processed chunk {i + 1}, found {len(chunk_results)} results"
                )
            else:
                print(f"  No results from chunk {i + 1}")

        except Exception as e:
            print(f"  Error processing chunk {i + 1}: {e}")
            continue

    # Combine all results
    if all_results:
        print(f"\nCombining results from {len(all_results)} chunks...")
        combined_results = pd.concat(all_results, ignore_index=True)

        # Add results to original mdata
        mdata.uns["test_results"] = combined_results
        print(f"Combined results: {len(combined_results)} total gene-element pairs")
    else:
        print("No results to combine!")
        mdata.uns["test_results"] = pd.DataFrame()

    # Save the updated mdata with covariates and results
    print(f"Saving results to {mdata_output_fp}...")
    mdata.write(mdata_output_fp, compression="gzip")

    return mdata


def run_perturbo_on_chunk(
    chunk_mdata,
    efficiency_mode,
    fit_guide_efficacy,
    accelerator,
    batch_size,
    early_stopping,
    early_stopping_patience,
    lr,
    num_epochs,
    gene_modality_name,
    guide_modality_name,
    test_all_pairs,
    control_guides,
):
    """
    Run PerTurbo inference on a single chunk.

    Returns:
    --------
    pd.DataFrame or None
        Results dataframe for this chunk, or None if no results
    """

    try:
        # Setup MuData for PerTurbo
        perturbo.PERTURBO.setup_mudata(
            chunk_mdata,
            perturbation_layer="guide_assignment",
            batch_key="batch",
            library_size_key="total_gene_umis",
            continuous_covariates_keys=["log1p_total_guide_umis_centered"],
            guide_by_element_key="intended_targets",
            gene_by_element_key="intended_targets" if not test_all_pairs else None,
            modalities={
                "rna_layer": gene_modality_name,
                "perturbation_layer": guide_modality_name,
            },
        )

        model = perturbo.PERTURBO(
            chunk_mdata,
            # control_guides=control_guides,  # broken in current PerTurbo version
            likelihood="nb",
            efficiency_mode=efficiency_mode,
            fit_guide_efficacy=fit_guide_efficacy,
        )

        model.train(
            num_epochs,
            lr=lr,
            batch_size=batch_size,
            accelerator=accelerator,
            early_stopping=early_stopping,
            early_stopping_patience=early_stopping_patience,
            early_stopping_min_delta=1e-5,
            early_stopping_monitor="elbo_train",
        )

        # Get element effects
        igvf_name_map = {
            "element": "intended_target_name",
            "gene": "gene_id",
            "q_value": "p_value",
        }

        element_effects = (
            model.get_element_effects()
            .rename(columns=igvf_name_map)
            .assign(log2_fc=lambda x: x["loc"] / np.log(2))
            .assign(
                gene_id=lambda x: x["gene_id"].astype("category"),
                intended_target_name=lambda x: x["intended_target_name"].astype(
                    "category"
                ),
            )
        )

        results = element_effects[
            [
                "gene_id",
                "intended_target_name",
                "log2_fc",
                "p_value",
            ]
        ]

        # Merge with pairs_to_test to get guide_id information
        if not test_all_pairs and "pairs_to_test" in chunk_mdata.uns:
            pairs_to_test_df = chunk_mdata.uns["pairs_to_test"]
            results = results.merge(
                pairs_to_test_df,
                on=["gene_id", "intended_target_name"],
                how="left",
            )
        else:
            # For test_all_pairs case
            results = results.merge(
                chunk_mdata[guide_modality_name].var[
                    ["intended_target_name", "guide_id"]
                ],
                how="left",
                on=["intended_target_name"],
            )

        # Rename columns to match expected output
        results.rename(
            columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"},
            inplace=True,
        )

        return results

    except Exception as e:
        print(f"    Error in PerTurbo processing: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Run chunked PerTurbo analysis on MuData",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument("mdata_output_fp", type=str, help="Output file path for MuData")

    parser.add_argument(
        "--chunk-size", "-c", type=int, default=4000, help="Number of genes per chunk"
    )

    parser.add_argument(
        "--fit_guide_efficacy",
        type=bool,
        default=True,
        help="Whether to fit guide efficacy (overrides efficiency_mode if false)",
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
        default="auto",
        help="Accelerator to use for training: default 'auto' (detects GPU if available on Linux, otherwise CPU)",
    )

    parser.add_argument(
        "--batch_size",
        type=int,
        default=512,
        help="Batch size for training",
    )

    parser.add_argument(
        "--early_stopping",
        type=bool,
        default=True,
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
        "--test_all_pairs",
        type=bool,
        default=True,
        help="Whether to test all pairs or only those in pairs_to_test",
    )

    parser.add_argument(
        "--test_control_guides",
        type=bool,
        default=True,
        help="Whether to remove control guides from the analysis",
    )

    parser.add_argument(
        "--num_workers",
        type=int,
        default=0,
        help="Number of workers for data loading",
    )

    args = parser.parse_args()

    # Validate arguments
    if not os.path.exists(args.mdata_input_fp):
        parser.error(f"Input file does not exist: {args.mdata_input_fp}")

    try:
        run_perturbo_chunked(
            mdata_input_fp=args.mdata_input_fp,
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
            test_all_pairs=args.test_all_pairs,
            test_control_guides=args.test_control_guides,
            num_workers=args.num_workers,
        )

        # print("Chunked PerTurbo inference completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
