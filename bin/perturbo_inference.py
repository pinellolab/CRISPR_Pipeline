#!/usr/bin/env python
import argparse
import perturbo
import mudata as md
import numpy as np
import pandas as pd


def run_perturbo(
    mdata_input_fp,
    mdata_output_fp,
    fit_guide_efficacy=True,  # whether to fit guide efficacy (if false, overrides efficiency_mode)
    efficiency_mode="auto",  # can be "mixture" (for low MOI only), "scaled", or "auto" (mixture for low MOI, scaled for high MOI)
    accelerator="gpu",  # can be "auto", "gpu" or "cpu"
    batch_size=512,  # batch size for training
    early_stopping=True,  # whether to use early stopping
    early_stopping_patience=5,  # patience for early stopping
    lr=0.01,  # learning rate for training
    num_epochs=100,  # (max) number of epochs for training
    gene_modality_name="gene",  # name of the gene modality in the MuData object
    guide_modality_name="guide",  # name of the guide modality in the MuData
):
    mdata = md.read(mdata_input_fp)
    mdata[gene_modality_name].obs = (
        mdata.obs.join(
            mdata[gene_modality_name].obs.drop(
                columns=mdata.obs.columns, errors="ignore"
            )
        )
        .join(
            mdata[guide_modality_name].obs.drop(
                columns=mdata.obs.columns.union(mdata[gene_modality_name].obs.columns),
                errors="ignore",
            )
        )
        .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"]))
    )

    mdata[guide_modality_name].X = mdata[guide_modality_name].layers["guide_assignment"]

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

    pairs_to_test_df = pd.DataFrame(mdata.uns["pairs_to_test"])
    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )
    # direct aggregate
    aggregated_df = (
        pairs_to_test_df.assign(value=1)
        .groupby(["gene_id", "intended_target_name"])
        .agg(value=("value", "max"))
        .reset_index()
    )

    # pivot the data
    mdata[gene_modality_name].varm["intended_targets"] = (
        aggregated_df.pivot(
            index="gene_id", columns="intended_target_name", values="value"
        )
        .reindex(mdata[gene_modality_name].var_names)
        .fillna(0)
    )

    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )

    intended_targets_df = pd.get_dummies(
        mdata[guide_modality_name].var["intended_target_name"]
    ).astype(float)

    mdata[guide_modality_name].varm["intended_targets"] = intended_targets_df[
        mdata.uns["intended_target_names"]
    ]

    ########################################
    ## subset mdata for perturbo speedup ##
    tested_guides = pairs_to_test_df["guide_id"].unique()
    tested_genes = pairs_to_test_df["gene_id"].unique()

    rna_subset = mdata[gene_modality_name][:, tested_genes]
    grna_feature_ids = mdata[guide_modality_name].var["guide_id"].isin(tested_guides)
    grna_subset = mdata[guide_modality_name][:, grna_feature_ids]

    mdata_dict = {gene_modality_name: rna_subset, guide_modality_name: grna_subset}
    if "hashing" in mdata.mod.keys():
        mdata_dict["hashing"] = mdata["hashing"]
    mdata_subset = md.MuData(mdata_dict)

    mdata_subset = mdata_subset.copy()
    ########################################

    perturbo.PERTURBO.setup_mudata(
        mdata_subset,
        batch_key="batch",
        library_size_key="total_gene_umis",
        continuous_covariates_keys=["log1p_total_guide_umis"],
        guide_by_element_key="intended_targets",
        gene_by_element_key="intended_targets",
        modalities={
            "rna_layer": gene_modality_name,
            "perturbation_layer": guide_modality_name,
        },
    )

    model = perturbo.PERTURBO(
        mdata_subset,
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

    igvf_name_map = {
        "element": "intended_target_name",
        "gene": "gene_id",
        "q_value": "p_value",
    }

    element_effects = (
        model.get_element_effects()
        .rename(columns=igvf_name_map)
        .assign(log2_fc=lambda x: x["loc"] / np.log(2))
        .merge(pairs_to_test_df)
    )

    mdata = md.read(mdata_input_fp)
    mdata.uns["test_results"] = element_effects[
        [
            "gene_id",
            "guide_id",
            "intended_target_name",
            "log2_fc",
            "p_value",
            "pair_type",
        ]
    ]

    mdata.uns["test_results"].rename(
        columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"},
        inplace=True,
    )

    mdata.write(mdata_output_fp)
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
        choices=["mixture", "scaled", "auto"],
        default="auto",
        help="Efficiency mode for the model (default: auto)",
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
        default=512,
        help="Batch size for training (default: 512)",
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
    )


if __name__ == "__main__":
    main()
