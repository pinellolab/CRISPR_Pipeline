#!/usr/bin/env python
import argparse
import perturbo
import mudata as md
import numpy as np
import pandas as pd
import scvi
from intended_target_key_utils import (
    annotate_intended_target_groups,
    enrich_pairs_with_target_metadata,
    get_target_lookup,
)


def _read_pairs_table(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith(".parquet"):
        return pd.read_parquet(path)
    if lower.endswith(".tsv") or lower.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    if lower.endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def run_perturbo(
    mdata_input_fp,
    results_tsv_fp,
    mdata_output_fp=None,
    fit_guide_efficacy=True,  # whether to fit guide efficacy (if false, overrides efficiency_mode)
    efficiency_mode="undecided",  # mapping from undecided->auto, low->mixture, high->scaled# can be "mixture" (for low MOI only), "scaled", "undecided" (auto), "low" (mixture), or "high" (scaled)
    accelerator="gpu",  # can be "auto", "gpu" or "cpu"
    batch_size=None,  # batch size for training
    early_stopping=False,  # whether to use early stopping
    early_stopping_patience=5,  # patience for early stopping
    lr=0.01,  # learning rate for training
    num_epochs=100,  # (max) number of epochs for training
    gene_modality_name="gene",  # name of the gene modality in the MuData object
    guide_modality_name="guide",  # name of the guide modality in the MuData
    test_all_pairs=False,  # whether to test all pairs or only those in pairs_to_test
    inference_type="element",  # can be per-guide or per-element
    drop_ntc_guides=False,  # whether to drop non-targeting control guides in low_MOI analysis
    num_workers=0,  # number of worker processes for data loading
    pairs_to_test_file=None,  # optional sidecar file with pairs_to_test
):
    scvi.settings.seed = 0
    if num_workers > 0:
        scvi.settings.dl_num_workers = num_workers

    mdata = md.read(mdata_input_fp)
    guide_var = annotate_intended_target_groups(mdata["guide"].var)
    mdata["guide"].var = guide_var

    if inference_type == "guide":
        element_key = "guide_id"
    elif inference_type == "element":
        element_key = "intended_target_key"
    else:
        raise ValueError("inference_type must be 'guide' or 'element'")

    guide_var = mdata["guide"].var
    target_lookup = get_target_lookup(guide_var) if inference_type == "element" else None

    control_guide_filter = pd.Series(False, index=guide_var.index)
    if "targeting" in guide_var.columns:
        targeting_series = guide_var["targeting"]
        if targeting_series.dtype != bool:
            targeting_series = (
                targeting_series.astype(str)
                .str.lower()
                .isin(["true", "1", "t", "yes", "y"])
            )
        control_guide_filter |= ~targeting_series
    if "type" in guide_var.columns:
        control_guide_filter |= (
            guide_var["type"].astype(str).str.lower().eq("non-targeting")
        )
    if "intended_target_name" in guide_var.columns:
        control_guide_filter |= (
            guide_var["intended_target_name"]
            .astype(str)
            .str.contains("non-targeting", case=False, na=False)
        )

    if not any(
        [c in guide_var.columns for c in ["targeting", "type", "intended_target_name"]]
    ):
        raise KeyError(
            "guide.var is missing all of: 'targeting', 'type', 'intended_target_name'. "
            "Cannot identify control guides."
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

    # efficiency_mode = {"undecided": "auto", "low": "mixture", "high": "scaled"}[
    #     efficiency_mode
    # ]
    efficiency_mode = "scaled"
    fit_guide_efficacy = True
    # if efficiency_mode == "auto":
    guides_per_element = mdata[guide_modality_name].var[element_key].value_counts()

    # max_guides_per_cell = mdata[guide_modality_name].X.sum(axis=1).max()
    if np.all(guides_per_element <= 1) or inference_type == "guide":
        fit_guide_efficacy = False
        print("Not fitting guide efficiency -- only one guide per element.")

    intended_targets_df = pd.get_dummies(
        mdata[guide_modality_name].var[element_key]
    ).astype(float)

    if efficiency_mode == "mixture":
        raise NotImplementedError(
            "Mixture efficiency mode is not currently supported due to issues with model convergence. Please use 'scaled' or 'undecided' instead."
        )
        # don't test for control guides in low_MOI analysis (slightly more robust alternative to just dropping "non-targeting")
        if drop_ntc_guides:
            control_elements_idx = (
                intended_targets_df.loc[control_guides].sum(axis=0) > 0
            )
            control_elements = intended_targets_df.columns[
                control_elements_idx
            ].tolist()
            intended_targets_df = intended_targets_df.drop(control_elements, axis=1)

        # filter any cells with >1 guide
        multi_guide_cells = (
            mdata[guide_modality_name].layers["guide_assignment"].sum(axis=1) > 1
        )
        if multi_guide_cells.any():
            print(
                f"Removing {multi_guide_cells.sum()} cells with multiple guides. ({multi_guide_cells.sum() / len(mdata) * 100:.1f}% of total)"
            )
            mdata = mdata[~multi_guide_cells, :].copy()

    mdata[guide_modality_name].varm["intended_targets"] = intended_targets_df
    mdata.uns[element_key] = intended_targets_df.columns.tolist()

    # create element by gene matrix if not testing all pairs
    if not test_all_pairs:
        if pairs_to_test_file is not None:
            pairs_to_test_df = _read_pairs_table(pairs_to_test_file)
        elif "pairs_to_test" in mdata.uns and isinstance(mdata.uns["pairs_to_test"], pd.DataFrame):
            pairs_to_test_df = mdata.uns["pairs_to_test"]
        elif "pairs_to_test" in mdata.uns and isinstance(mdata.uns["pairs_to_test"], dict):
            pairs_to_test_df = pd.DataFrame(mdata.uns["pairs_to_test"])
        else:
            raise ValueError(
                "pairs_to_test not found (provide --pairs_to_test_file or include pairs_to_test in mudata.uns)."
            )

        if "gene_id" not in pairs_to_test_df.columns and "gene_name" in pairs_to_test_df.columns:
            pairs_to_test_df = pairs_to_test_df.rename(columns={"gene_name": "gene_id"})
        if "gene_id" not in pairs_to_test_df.columns:
            raise ValueError("pairs_to_test must contain gene_id (or gene_name) column.")

        if element_key not in pairs_to_test_df.columns:
            pairs_to_test_df = enrich_pairs_with_target_metadata(pairs_to_test_df, guide_var)

        aggregated_df = (
            pairs_to_test_df[["gene_id", element_key]].drop_duplicates().assign(value=1)
        )

        # pivot the data
        mdata[gene_modality_name].varm["intended_targets"] = (
            aggregated_df.pivot(index="gene_id", columns=element_key, values="value")
            .reindex(
                index=mdata[gene_modality_name].var_names,
                columns=mdata.uns[element_key],
            )
            .fillna(0)
        )

    ########################################
    # Setup MuData for PerTurbo

    perturbo.PERTURBO.setup_mudata(
        mdata,
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

    if control_guides is not None and isinstance(control_guides[0], str):
        control_guides = mdata["guide"].var_names.isin(control_guides)

    model = perturbo.PERTURBO(
        mdata,
        control_guides=control_guides,
        likelihood="nb",
        efficiency_mode=efficiency_mode,
        fit_guide_efficacy=fit_guide_efficacy,
    )

    model.view_anndata_setup(mdata, hide_state_registries=True)
    if batch_size is None:
        batch_size = int(np.clip(len(mdata) // 20, 512, 10_000))

    print(f"Training using batch size of {batch_size}")
    model.train(
        num_epochs,  # max number of epochs
        lr=lr,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_min_delta=1e-5,
        early_stopping_monitor="elbo_train",
        # data_splitter_kwargs={"drop_last": True},  # requires scvi-tools version > 1.2
    )

    # Reformat the output to match IGVF specifications
    igvf_name_map = {
        "element": element_key,
        "gene": "gene_id",
        "q_value": "p_value",
    }

    element_effects = (
        model.get_element_effects()
        .rename(columns=igvf_name_map)
        .assign(log2_fc=lambda x: x["loc"] / np.log(2))
        .assign(log2_scale=lambda x: x["scale"] / np.log(2))
    )

    element_effects.to_csv(
        "perturbo_results.tsv.gz", index=False, sep="\t", compression="gzip"
    )

    # element_effects[element_key] = element_effects[element_key].astype("category")
    # element_effects["gene_id"] = element_effects["gene_id"].astype("category")

    test_results = element_effects[
        [
            "gene_id",
            element_key,
            "log2_fc",
            "p_value",
        ]
    ]

    if inference_type == "element":
        test_results = test_results.merge(
            target_lookup,
            on="intended_target_key",
            how="left",
        )
        if test_results["intended_target_name"].isna().any():
            missing_keys = (
                test_results.loc[
                    test_results["intended_target_name"].isna(), "intended_target_key"
                ]
                .astype(str)
                .drop_duplicates()
                .tolist()
            )
            raise ValueError(
                "Unable to map intended_target_key back to target metadata for keys: "
                + ", ".join(missing_keys[:20])
            )
        test_results = test_results[
            [
                "gene_id",
                "intended_target_name",
                "intended_target_chr",
                "intended_target_start",
                "intended_target_end",
                "log2_fc",
                "p_value",
            ]
        ]

    mdata.uns[f"per_{inference_type}_results"] = test_results

    # NOTE: this part creates a per-guide output table even when we are running per-element inference.
    # This is to maintain compatibility with the existing workflow, which condenses per-guide output
    # into per-element output in a separate module.

    # if not test_all_pairs:
    #     mdata.uns["test_results"] = mdata.uns["test_results"].merge(
    #         pairs_to_test_df,
    #         on=["gene_id", element_key],
    #         how="left",
    #     )
    # else:
    #     mdata.uns["test_results"] = mdata.uns["test_results"].merge(
    #         mdata["guide"].var[["intended_target_name", "guide_id"]],
    #         how="left",
    #         on=[element_key],
    #     )

    # mdata.uns["test_results"].rename(
    #     columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"},
    #     inplace=True,
    # )

    # Write results table to TSV (required)
    print("Writing results to ", results_tsv_fp)
    if results_tsv_fp.endswith(".gz"):
        test_results.to_csv(results_tsv_fp, index=False, sep="\t", compression="gzip")
    else:
        test_results.to_csv(results_tsv_fp, index=False, sep="\t")

    # Optionally write the full MuData if an output path was provided
    if mdata_output_fp:
        print("Writing mudata to ", mdata_output_fp)
        mdata.write(mdata_output_fp, compression="gzip")

    return mdata


def main():
    parser = argparse.ArgumentParser(description="Run PerTurbo analysis on MuData")
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument(
        "results_tsv_fp",
        type=str,
        help="Output TSV file path for mdata.uns['test_results'] (required)",
    )
    parser.add_argument(
        "--mdata_output_fp",
        type=str,
        default=None,
        help="Optional output file path for MuData; if omitted the MuData will not be written",
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
        help="Efficiency mode for the model: 'undecided'/'auto' (auto-detect), 'low' (mixture), 'high' (scaled), 'mixture', or 'scaled' (default: undecided)",
    )

    parser.add_argument(
        "--accelerator",
        type=str,
        choices=["auto", "gpu", "cpu"],
        default="gpu",
        help="Accelerator to use for training (default: gpu)",
    )
    # parser.add_argument(
    #     "--batch_size",
    #     type=int,
    #     default=4096,
    #     help="Batch size for training (default: 4096)",
    # )
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
        "--inference_type",
        type=str,
        default="element",
        help="Unit to test for effects on each gene: 'guide' or 'element' (default: 'element')",
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
    parser.add_argument(
        "--pairs_to_test_file",
        type=str,
        default=None,
        help="Optional sidecar pairs table path (tsv/csv/parquet) for non-all-pairs runs.",
    )

    # Parse the arguments
    args = parser.parse_args()

    run_perturbo(
        args.mdata_input_fp,
        args.results_tsv_fp,
        mdata_output_fp=args.mdata_output_fp,
        fit_guide_efficacy=args.fit_guide_efficacy,
        efficiency_mode=args.efficiency_mode,
        accelerator=args.accelerator,
        # batch_size=args.batch_size,
        early_stopping=args.early_stopping,
        early_stopping_patience=args.early_stopping_patience,
        lr=args.lr,
        num_epochs=args.num_epochs,
        gene_modality_name=args.gene_modality_name,
        guide_modality_name=args.guide_modality_name,
        test_all_pairs=args.test_all_pairs,
        inference_type=args.inference_type,
        pairs_to_test_file=args.pairs_to_test_file,
    )


if __name__ == "__main__":
    main()
