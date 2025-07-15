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
    efficiency_mode="mixture",  # can be "mixture" (for low MOI) or "scaled" (low or high MOI)
    accelerator="gpu",  # can be "auto", "gpu" or "cpu"
    batch_size=512,  # batch size for training
    early_stopping=True,  # whether to use early stopping
    early_stopping_patience=5,  # patience for early stopping
    lr=0.01,  # learning rate for training
    num_epochs=100,  # (max) number of epochs for training
):
    mdata = md.read(mdata_input_fp)
    mdata["gene"].obs = (
        mdata.obs.join(
            mdata["gene"].obs.drop(columns=mdata.obs.columns, errors="ignore")
        )
        .join(
            mdata["guide"].obs.drop(
                columns=mdata.obs.columns.union(mdata["gene"].obs.columns),
                errors="ignore",
            )
        )
        .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"]))
    )

    mdata["guide"].X = mdata["guide"].layers["guide_assignment"]
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
    mdata["gene"].varm["intended_targets"] = (
        aggregated_df.pivot(
            index="gene_id", columns="intended_target_name", values="value"
        )
        .reindex(mdata["gene"].var_names)
        .fillna(0)
    )

    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )

    intended_targets_df = pd.get_dummies(
        mdata["guide"].var["intended_target_name"]
    ).astype(float)

    mdata["guide"].varm["intended_targets"] = intended_targets_df[
        mdata.uns["intended_target_names"]
    ]

    ########################################
    ## subset mdata for perturbo speedup ##
    tested_guides = pairs_to_test_df["guide_id"].unique()
    tested_genes = pairs_to_test_df["gene_id"].unique()

    rna_subset = mdata["gene"][:, tested_genes]
    grna_feature_ids = mdata["guide"].var["guide_id"].isin(tested_guides)
    grna_subset = mdata["guide"][:, grna_feature_ids]

    mdata_dict = {"gene": rna_subset, "guide": grna_subset}
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
            "rna_layer": "gene",
            "perturbation_layer": "guide",
        },
    )

    model = perturbo.PERTURBO(
        mdata_subset,
        likelihood="nb",
        efficiency_mode=efficiency_mode,
    )

    model.train(
        num_epochs=num_epochs,
        lr=lr,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        fit_guide_efficacy=fit_guide_efficacy,
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
        action="store_true",
        help="Whether to fit guide efficacy (overrides efficiency_mode if false)",
    )
    parser.add_argument(
        "--efficiency_mode",
        type=str,
        choices=["mixture", "scaled"],
        default="scaled",
        help="Efficiency mode for the model (default: scaled)",
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
        action="store_true",
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
    )


if __name__ == "__main__":
    main()
