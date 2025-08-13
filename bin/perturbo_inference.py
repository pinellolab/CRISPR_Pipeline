#!/usr/bin/env python
import argparse
import perturbo
import mudata as md
import numpy as np
import pandas as pd
import scvi


def run_perturbo(
    mdata_input_fp,
    mdata_output_fp,
    fit_guide_efficacy=True,  # whether to fit guide efficacy (if false, overrides efficiency_mode)
    efficiency_mode="undecided",  # mapping from undecided->auto, low->mixture, high->scaled# can be "mixture" (for low MOI only), "scaled", "undecided" (auto), "low" (mixture), or "high" (scaled)
    accelerator="gpu",  # can be "auto", "gpu" or "cpu"
    batch_size=4096,  # batch size for training
    early_stopping=False,  # whether to use early stopping
    early_stopping_patience=5,  # patience for early stopping
    lr=0.01,  # learning rate for training
    num_epochs=100,  # (max) number of epochs for training
    gene_modality_name="gene",  # name of the gene modality in the MuData object
    guide_modality_name="guide",  # name of the guide modality in the MuData
    test_all_pairs=False,  # whether to test all pairs or only those in pairs_to_test
    test_control_guides=True,  # whether to remove control guides from the analysis
    num_workers=0,  # number of worker processes for data loading
):
    scvi.settings.seed = 0
    if num_workers > 0:
        scvi.settings.dl_num_workers = num_workers

    mdata = md.read(mdata_input_fp)

    control_guide_filter = (~mdata["guide"].var["targeting"]) | (
        mdata["guide"]
        .var["type"]
        .isin(["safe-targeting", "non-targeting", "negative control"])
    )
    if np.any(control_guide_filter):
        control_guides = mdata["guide"].var_names[control_guide_filter].tolist()
    else:
        control_guides = None

    mdata[gene_modality_name].obs["log1p_total_guide_umis"] = np.log1p(
        mdata[guide_modality_name].obs["total_guide_umis"]
    )
    mdata[gene_modality_name].obs["log1p_total_guide_umis_centered"] = (
        mdata[gene_modality_name].obs["log1p_total_guide_umis"]
        - mdata[gene_modality_name].obs["log1p_total_guide_umis"].mean()
    )

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

    intended_targets_df = pd.get_dummies(
        mdata[guide_modality_name].var["intended_target_name"]
    ).astype(float)
    if not test_control_guides:
        control_elements = intended_targets_df[control_guides].sum(axis=0) > 0
        intended_targets_df = intended_targets_df.drop(control_elements, axis=1)

    mdata[guide_modality_name].varm["intended_targets"] = intended_targets_df

    mdata.uns["intended_target_names"] = intended_targets_df.columns.tolist()

    # create element by gene matrix if not testing all pairs
    if not test_all_pairs:
        if isinstance(mdata.uns["pairs_to_test"], pd.DataFrame):
            pairs_to_test_df = mdata.uns["pairs_to_test"]
        elif isinstance(mdata.uns["pairs_to_test"], dict):
            pairs_to_test_df = pd.DataFrame(mdata.uns["pairs_to_test"])

        aggregated_df = (
            pairs_to_test_df[["gene_id", "intended_target_name"]]
            .drop_duplicates()
            .assign(value=1)
        )

        # pivot the data
        mdata[gene_modality_name].varm["intended_targets"] = (
            aggregated_df.pivot(
                index="gene_id", columns="intended_target_name", values="value"
            )
            .reindex(
                index=mdata[gene_modality_name].var_names,
                columns=mdata.uns["intended_target_names"],
            )
            .fillna(0)
        )

        # subset mdata for perturbo speedup
        tested_guides = pairs_to_test_df["guide_id"].unique()
        tested_genes = pairs_to_test_df["gene_id"].unique()

        rna_subset = mdata[gene_modality_name][:, tested_genes]
        grna_feature_ids = (
            mdata[guide_modality_name].var["guide_id"].isin(tested_guides)
        )
        grna_subset = mdata[guide_modality_name][:, grna_feature_ids]

        targeted_cells = grna_subset.X.sum(axis=1) > 0
        if not targeted_cells.any():
            raise ValueError("No targeted cells found in the guide modality subset.")

        mdata_dict = {
            gene_modality_name: rna_subset[targeted_cells],
            guide_modality_name: grna_subset[targeted_cells],
        }

        # copy over any additional modalities
        for mod in mdata.mod.keys():
            if mod not in mdata_dict:
                mdata_dict[mod] = mdata[mod][targeted_cells]
        mdata_subset = md.MuData(mdata_dict).copy()
    else:
        mdata_subset = mdata

    ########################################
    # Setup MuData for PerTurbo

    perturbo.PERTURBO.setup_mudata(
        mdata_subset,
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
        mdata_subset,
        # control_guides=control_guides, # broken in current PerTurbo version, fix when we update image
        likelihood="nb",
        efficiency_mode=efficiency_mode,
        fit_guide_efficacy=fit_guide_efficacy,
    )

    model.train(
        num_epochs,  # max number of epochs
        lr=lr,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_min_delta=1e-5,
        early_stopping_monitor="elbo_train",
    )

    # Reformat the output to match IGVF specifications
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
            intended_target_name=lambda x: x["intended_target_name"].astype("category"),
        )
    )

    mdata.uns["test_results"] = element_effects[
        [
            "gene_id",
            "intended_target_name",
            "log2_fc",
            "p_value",
        ]
    ]

    # NOTE: this part creates a per-guide output table even though we are running per-element inference.
    # This is to maintain compatibility with the existing workflow, which condenses per-guide output
    # into per-element output in a separate module.

    if not test_all_pairs:
        mdata.uns["test_results"] = mdata.uns["test_results"].merge(
            pairs_to_test_df,
            on=["gene_id", "intended_target_name"],
            how="left",
        )
    else:
        mdata.uns["test_results"] = mdata.uns["test_results"].merge(
            mdata["guide"].var[["intended_target_name", "guide_id"]],
            how="left",
            on=["intended_target_name"],
        )

    mdata.uns["test_results"].rename(
        columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"},
        inplace=True,
    )

    mdata.write(mdata_output_fp, compression="gzip")
    return mdata


def main():
    parser = argparse.ArgumentParser(description="Run PerTurbo analysis on MuData")
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument("mdata_output_fp", type=str, help="Output file path for MuData")
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
        help="Efficiency mode for the model: 'undecided'/'auto' (auto-detect), 'low' (mixture), 'high' (scaled), 'mixture', or 'scaled' (default: undecided)",
    )

    parser.add_argument(
        "--accelerator",
        type=str,
        choices=["auto", "gpu", "cpu"],
        default="gpu",
        help="Accelerator to use for training (default: gpu)",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=4096,
        help="Batch size for training (default: 4096)",
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
        help="Patience for early stopping (default: 5)",
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=0.01,
        help="Learning rate for training (default: 0.01)",
    )
    parser.add_argument(
        "--num_epochs",
        type=int,
        default=100,
        help="Maximum number of epochs for training (default: 100)",
    )
    parser.add_argument(
        "--gene_modality_name",
        type=str,
        default="gene",
        help="Name of the gene modality in the MuData object (default: 'gene')",
    )
    parser.add_argument(
        "--guide_modality_name",
        type=str,
        default="guide",
        help="Name of the guide modality in the MuData object (default: 'guide')",
    )
    parser.add_argument(
        "--test_all_pairs",
        action="store_true",
        help="Whether to test all pairs or only those in pairs_to_test (default: False)",
    )
    parser.add_argument(
        "--test_control_guides",
        type=bool,
        default=True,
        help="Whether to remove control guides from the analysis (default: False)",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=0,
        help="Number of workers for data loading (default: 0)",
    )

    # Parse the arguments
    args = parser.parse_args()
    run_perturbo(
        args.mdata_input_fp,
        args.mdata_output_fp,
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
    )


if __name__ == "__main__":
    main()
