#!/usr/bin/env python
"""Run PerTurbo v2 while preserving CRISPR_Pipeline's inference outputs.

The score-based conditional randomization test builds on SCEPTRE (Barry et al.,
2024, https://doi.org/10.1186/s13059-024-03254-2) and score-resampling work
(Barry et al., 2025, https://arxiv.org/abs/2501.03530). Its Bernoulli saddlepoint
approximation builds on spaCRT (Niu et al., https://arxiv.org/abs/2407.08911).
PerTurbo provides the GPU implementation and integration with Bayesian effect
estimation.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import gzip

import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control

from analysis_output_formatting import make_h5mu_safe_dataframe
from control_group import SETTING_PARAM, normalize_moi
from intended_target_key_utils import (
    annotate_intended_target_groups,
    enrich_pairs_with_target_metadata,
    get_target_lookup,
)
from mudata_uns_io import write_uns_patch
from inference_covariates import (
    CANONICAL_COVARIATES,
    GUIDE_UMI_COLUMN,
    derive_perturbo_log_covariates,
    ensure_gene_depth_totals,
    ensure_guide_umi_totals,
    perturbo_covariate_arguments,
    perturbo_log_name,
)

GUIDE_UMI_COVARIATE = perturbo_log_name(GUIDE_UMI_COLUMN)
from result_table_io import write_result_table
from bounded_perturbo_results import parquet_manifest, sha256_file, verify_manifest, write_bh_sidecar, write_manifest
from compact_result_dtypes import compact_result_floats


GENE_MODALITY = "gene"
GUIDE_MODALITY = "guide"
ASSIGNMENT_LAYER = "guide_assignment"
ELEMENT_MAP_KEY = "perturbo_v2_intended_targets"
ELEMENT_NAMES_KEY = "perturbo_v2_intended_target_names"
GUIDE_MAP_KEY = "perturbo_v2_guides"
GUIDE_NAMES_KEY = "perturbo_v2_guide_names"
CONTROL_SUBSTRING = "non-targeting"


@contextmanager
def _open_mudata(path: str | Path, *, backed: str | None = None):
    """Open MuData and deterministically release its HDF5 file handle."""
    mdata = md.read_h5mu(path, backed=backed)
    try:
        yield mdata
    finally:
        file_manager = getattr(mdata, "file", None)
        if file_manager is not None:
            file_manager.close()


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if not valid.any():
        return out
    p_valid = p.loc[valid].clip(lower=0.0, upper=1.0)
    out.loc[p.loc[valid].index] = false_discovery_control(
        p_valid.to_numpy(dtype=float), method="bh"
    )
    return out


def _guide_control_mask(guide_var: pd.DataFrame) -> pd.Series:
    targeting = guide_var.get("targeting")
    if targeting is None:
        targeting_mask = pd.Series(False, index=guide_var.index)
    elif targeting.dtype == bool:
        targeting_mask = targeting.fillna(False)
    else:
        targeting_mask = (
            targeting.astype(str).str.strip().str.lower().isin({"true", "1", "t", "yes", "y"})
        )

    guide_type = guide_var.get("type", pd.Series("", index=guide_var.index))
    guide_type = guide_type.astype(str).str.strip().str.lower()
    target_name = guide_var.get("intended_target_name", pd.Series("", index=guide_var.index))
    target_name = target_name.astype(str).str.strip().str.lower()
    return (
        ~targeting_mask
        | guide_type.isin({"safe-targeting", "non-targeting", "negative control"})
        | target_name.eq("non-targeting")
        | target_name.str.startswith("non-targeting|")
    )


def _get_assignment_matrix(mdata: md.MuData):
    guide = mdata[GUIDE_MODALITY]
    matrix = guide.layers[ASSIGNMENT_LAYER] if ASSIGNMENT_LAYER in guide.layers else guide.X
    return matrix if sparse.issparse(matrix) else sparse.csr_matrix(matrix)


def _ensure_covariates(mdata: md.MuData) -> None:
    """Take the PerTurbo log covariates from upstream, or derive them once here.

    The derivations live in ``bin/inference_covariates.py``; the shared
    preparation step already ran them over the analysed cells, so when the
    columns are present they are used as they stand. Recomputing here is only for
    input paths that skip that step -- and the underlying totals are pinned first
    so that a depth is never measured over an already-subset matrix.
    """
    expected = [
        c.perturbo_name for c in CANONICAL_COVARIATES if c.kind == "count"
    ]
    if all(name in mdata[GENE_MODALITY].obs.columns for name in expected):
        print(f"Using the upstream PerTurbo log covariates: {expected}.")
        return
    ensure_guide_umi_totals(mdata)
    ensure_gene_depth_totals(mdata)
    derive_perturbo_log_covariates(mdata)


def _ensure_library_size(mdata: md.MuData) -> None:
    gene = mdata[GENE_MODALITY]
    if "total_gene_umis" not in gene.obs.columns:
        gene.obs["total_gene_umis"] = np.asarray(gene.X.sum(axis=1)).ravel()


def _build_element_mapping(mdata: md.MuData) -> pd.DataFrame:
    guide = mdata[GUIDE_MODALITY]
    guide.var = annotate_intended_target_groups(guide.var)
    mapping = pd.get_dummies(guide.var["intended_target_key"]).astype(float)
    guide.varm[ELEMENT_MAP_KEY] = mapping
    guide.uns[ELEMENT_NAMES_KEY] = mapping.columns.astype(str).tolist()
    return mapping


def _build_guide_identity_mapping(mdata: md.MuData) -> dict[str, str]:
    guide = mdata[GUIDE_MODALITY]
    guide_var = guide.var.copy()
    if "guide_id" not in guide_var.columns:
        guide_var["guide_id"] = guide_var.index.astype(str)
    guide_ids = guide_var["guide_id"].astype(str)
    control_mask = _guide_control_mask(guide_var)
    perturbo_names = guide_ids.copy()
    perturbo_names.loc[control_mask] = CONTROL_SUBSTRING + "|" + guide_ids.loc[control_mask]

    mapping = sparse.eye(len(guide_ids), dtype=np.float32, format="csr")
    guide.varm[GUIDE_MAP_KEY] = mapping
    guide.uns[GUIDE_NAMES_KEY] = perturbo_names.astype(str).tolist()
    return dict(zip(perturbo_names.astype(str), guide_ids.tolist()))


def _covariate_has_baseline_variance(
    input_path: Path, map_key: str, names_key: str, *, pool: str
) -> bool:
    """Does the guide-UMI covariate vary among the cells the baseline is fit on?

    Which cells those are is the CRT pool's business, so the check has to follow
    it. Under ``control-anchored`` PerTurbo fits its stage-one baseline on the
    pure control cells, and a covariate constant across them is unidentifiable
    there. Under ``all-cells`` -- the shipped high-MOI default -- the baseline is
    fit on every cell, so the control cells say nothing about identifiability:
    measuring variance over the controls alone would drop the covariate from a
    screen that merely has few or no non-targeting guides, and PerTurbo would
    then be conditioning on less than SCEPTRE, which is the asymmetry this
    covariate list exists to remove. ``auto`` leaves the choice to PerTurbo, so
    it is treated as the restrictive case.
    """
    with _open_mudata(input_path, backed="r") as mdata:
        if GUIDE_UMI_COVARIATE not in mdata[GENE_MODALITY].obs.columns:
            return False
        values = pd.to_numeric(
            mdata[GENE_MODALITY].obs[GUIDE_UMI_COVARIATE],
            errors="coerce",
        ).to_numpy(dtype=float)
        if pool == "all-cells":
            finite = values[np.isfinite(values)]
            if finite.size < 2:
                return False
            return bool(np.nanstd(finite) > 1e-8)
        guide = mdata[GUIDE_MODALITY]
        mapping = guide.varm[map_key]
        names = [str(x) for x in guide.uns[names_key]]
        control_idx = [i for i, name in enumerate(names) if CONTROL_SUBSTRING in name.lower()]
        if not control_idx:
            return False
        if isinstance(mapping, pd.DataFrame):
            mapping = mapping.to_numpy()
        control_guides = np.asarray(mapping[:, control_idx].sum(axis=1)).ravel() > 0
        if not np.any(control_guides):
            return False
        assignment = _get_assignment_matrix(mdata)
        control_cells = np.asarray((assignment[:, control_guides] > 0).sum(axis=1)).ravel() > 0
        control_values = values[control_cells & np.isfinite(values)]
        if control_values.size < 2:
            return False
        return bool(np.nanstd(control_values) > 1e-8)


def _build_native_pairs_to_test(
    mdata: md.MuData,
    *,
    pair_element_column: str,
    element_names: list[str],
    guide_name_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build PerTurbo's native two-column ``element,gene`` restriction table."""
    pairs = mdata.uns.get("pairs_to_test")
    if pairs is None:
        raise KeyError("pairs_to_test not found in MuData; local PerTurbo fitting requires explicit pairs.")
    pairs = pairs.copy() if isinstance(pairs, pd.DataFrame) else pd.DataFrame(pairs)
    if pair_element_column == "intended_target_key" and pair_element_column not in pairs.columns:
        pairs = enrich_pairs_with_target_metadata(pairs, pd.DataFrame(mdata[GUIDE_MODALITY].var))

    gene_names = set(mdata[GENE_MODALITY].var_names.astype(str))
    element_name_set = set(element_names)
    if guide_name_map is None:
        pair_element_names = pairs[pair_element_column].astype(str)
    else:
        output_name_by_guide = {guide_id: output_name for output_name, guide_id in guide_name_map.items()}
        pair_element_names = pairs[pair_element_column].astype(str).map(output_name_by_guide)

    pair_gene_names = pairs["gene_id"].astype(str)
    missing_genes = sorted(set(pair_gene_names) - gene_names)
    missing_elements = sorted(set(pair_element_names.dropna()) - element_name_set)
    unmapped_elements = pairs.loc[pair_element_names.isna(), pair_element_column].astype(str).drop_duplicates().tolist()
    if missing_genes or missing_elements or unmapped_elements:
        messages = []
        if missing_genes:
            messages.append(f"{len(missing_genes)} genes absent from RNA modality: {', '.join(missing_genes[:10])}")
        if missing_elements:
            messages.append(f"{len(missing_elements)} elements absent from mapping: {', '.join(missing_elements[:10])}")
        if unmapped_elements:
            messages.append(f"{len(unmapped_elements)} guide IDs could not be mapped: {', '.join(unmapped_elements[:10])}")
        raise ValueError("Cannot construct PerTurbo local pairs-to-test table; " + "; ".join(messages))

    requested = pd.DataFrame(
        {
            "element": pair_element_names.to_numpy(dtype=str),
            "gene": pair_gene_names.to_numpy(dtype=str),
        }
    ).drop_duplicates(ignore_index=True)

    # The requested pairs alone, deliberately. Earlier versions added every
    # non-targeting element crossed with every tested gene, because the fit was
    # restricted to the requested pairs and PerTurbo calibrated empirical p-values
    # against fitted control-element z-values, so the controls had to be requested
    # too. Neither holds now: the fit covers every pair and the p-value comes from
    # the conditional randomization test. What the augmentation did do was put those
    # control pairs in this table's Benjamini-Hochberg family - on the Replogle
    # screen 2,623,521 of 2,714,942 rows, 96.6% - which made the local q-values far
    # more conservative than SCEPTRE's over the same hypotheses (SCEPTRE corrects
    # over the 91,421 requested pairs alone). Control pairs are still tested and
    # still present in the transcriptome-wide table, which is what the control
    # evaluation step reads.
    print(f"Prepared native PerTurbo pairs-to-test table: {len(requested):,} requested pairs.")
    return requested


def prepare_mudata_for_perturbo_v2(
    input_path: str | Path,
    output_path: str | Path,
    *,
    test_all_pairs: bool = False,
) -> dict[str, str]:
    """Write a PerTurbo v2-ready MuData and return guide output-name remapping."""
    with _open_mudata(input_path) as mdata:
        if GENE_MODALITY not in mdata.mod or GUIDE_MODALITY not in mdata.mod:
            raise KeyError("Expected MuData modalities named 'gene' and 'guide'.")
        _ensure_library_size(mdata)
        _ensure_covariates(mdata)
        _build_element_mapping(mdata)
        guide_name_map = _build_guide_identity_mapping(mdata)
        if not test_all_pairs and "pairs_to_test" not in mdata.uns:
            raise KeyError("pairs_to_test not found in MuData; local PerTurbo fitting requires explicit pairs.")
        mdata.write(output_path)
    return guide_name_map


def _maybe_filter_pairs(df: pd.DataFrame, prepared_mudata_path: str | Path, *, inference_type: str, test_all_pairs: bool) -> pd.DataFrame:
    if test_all_pairs:
        return df
    with _open_mudata(prepared_mudata_path, backed="r") as mdata:
        if "pairs_to_test" not in mdata.uns:
            raise KeyError("pairs_to_test not found in MuData; use --test-all-pairs to disable pair filtering.")
        pairs = mdata.uns["pairs_to_test"]
        if not isinstance(pairs, pd.DataFrame):
            pairs = pd.DataFrame(pairs)
        else:
            pairs = pairs.copy()
        guide_var = pd.DataFrame(mdata[GUIDE_MODALITY].var).copy()
    if inference_type == "element":
        if "intended_target_key" not in pairs.columns:
            pairs = enrich_pairs_with_target_metadata(pairs, guide_var)
        keep = pairs[["gene_id", "intended_target_key"]].drop_duplicates()
        return df.merge(keep, on=["gene_id", "intended_target_key"], how="inner")
    keep = pairs[["gene_id", "guide_id"]].drop_duplicates()
    return df.merge(keep, on=["gene_id", "guide_id"], how="inner")


def convert_element_effects(
    effects: pd.DataFrame,
    prepared_mudata_path: str | Path,
    *,
    test_all_pairs: bool,
    crt: bool = False,
) -> pd.DataFrame:
    with _open_mudata(prepared_mudata_path, backed="r") as mdata:
        target_lookup = get_target_lookup(pd.DataFrame(mdata[GUIDE_MODALITY].var).copy())
    out = _convert_common_effect_columns(effects, crt=crt).rename(columns={"element": "intended_target_key"})
    out = _maybe_filter_pairs(out, prepared_mudata_path, inference_type="element", test_all_pairs=test_all_pairs)
    out = out.merge(target_lookup, on="intended_target_key", how="left")
    missing = out["intended_target_name"].isna()
    if missing.any():
        keys = out.loc[missing, "intended_target_key"].astype(str).drop_duplicates().tolist()
        raise ValueError("Unable to map PerTurbo element names to target metadata: " + ", ".join(keys[:20]))
    out = out[
        [
            "gene_id",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
            "log2_fc",
            "perturbo_fc_se",
            "p_value",
            "perturbo_posterior_prob",
        ]
        + _carried_diagnostic_columns(out)
    ]
    out["perturbo_q_value"] = _bh_adjust(out["p_value"])
    return out


def convert_guide_effects(
    effects: pd.DataFrame,
    guide_name_map: dict[str, str],
    prepared_mudata_path: str | Path,
    *,
    test_all_pairs: bool,
    crt: bool = False,
) -> pd.DataFrame:
    out = _convert_common_effect_columns(effects, crt=crt)
    out["guide_id"] = out["element"].map(guide_name_map).fillna(out["element"]).astype(str)
    out = _maybe_filter_pairs(out, prepared_mudata_path, inference_type="guide", test_all_pairs=test_all_pairs)
    out = out[
        ["gene_id", "guide_id", "log2_fc", "perturbo_fc_se", "p_value", "perturbo_posterior_prob"]
        + _carried_diagnostic_columns(out)
    ]
    out["perturbo_q_value"] = _bh_adjust(out["p_value"])
    return out


# Per-pair CRT diagnostics PerTurbo writes beside the p-value. They travel
# into the pipeline's result tables under a `perturbo_` prefix so that a call
# can be read against how much data it rests on and how its tail probability
# was obtained, without opening PerTurbo's own artifact directory. None of them
# filters anything: `crt_low_information` is an annotation, not a gate, and the
# tail columns record why an approximation was replaced, not that the
# replacement is invalid. Whichever of them a PerTurbo version emits are carried;
# the tail-policy columns (`crt_tail_failure_reason` onward) arrive with the
# Chernoff-fallback release and are simply absent before it.
CRT_DIAGNOSTIC_COLUMNS: tuple[str, ...] = (
    "crt_low_information",
    "crt_observed_nonzero",
    "crt_expected_nonzero",
    "crt_saddlepoint_valid",
    "crt_tail_failure_reason",
    "crt_used_chernoff",
    "crt_used_conservative_one",
    "crt_root_residual_null_sd",
)


def _carried_diagnostic_columns(frame: pd.DataFrame) -> list[str]:
    return [f"perturbo_{c}" for c in CRT_DIAGNOSTIC_COLUMNS if f"perturbo_{c}" in frame.columns]


def _convert_common_effect_columns(effects: pd.DataFrame, *, crt: bool = False) -> pd.DataFrame:
    required = {"element", "gene", "posterior_mean", "posterior_scale", "posterior_prob"}
    missing = sorted(required - set(effects.columns))
    if missing:
        raise ValueError("PerTurbo element_effects.parquet is missing columns: " + ", ".join(missing))
    out = effects.copy()
    # The conditional randomization test's p-value is the primary significance
    # measure when the run produced one; the posterior probability is kept beside
    # it. Without the CRT the previous fallbacks apply unchanged.
    if "crt_saddlepoint_p_value" in out.columns:
        p_value = pd.to_numeric(out["crt_saddlepoint_p_value"], errors="coerce")
    elif "crt_p_value" in out.columns:
        p_value = pd.to_numeric(out["crt_p_value"], errors="coerce")
    elif "empirical_p_value" in out.columns:
        p_value = pd.to_numeric(out["empirical_p_value"], errors="coerce")
    else:
        p_value = pd.Series(np.nan, index=out.index, dtype=float)
    posterior = pd.to_numeric(out["posterior_prob"], errors="coerce")
    if crt:
        # No fallback when the conditional randomization test ran. It does not test
        # control elements - their cells are the pool it resamples within - so those
        # rows have no p-value, and filling them with the posterior probability put
        # two different quantities in one column: on the Replogle screen every one of
        # 2,623,521 control rows carried a posterior probability while the targeting
        # rows carried CRT p-values. Anything comparing the two, the pipeline's own
        # control evaluation included, was comparing incomparable numbers. A missing
        # p-value stays missing; Benjamini-Hochberg already preserves NaN.
        out["p_value"] = p_value
    else:
        out["p_value"] = p_value.where(p_value.notna(), posterior)
    out["perturbo_posterior_prob"] = posterior
    for column in CRT_DIAGNOSTIC_COLUMNS:
        if column in out.columns:
            out[f"perturbo_{column}"] = out[column]
    out["gene_id"] = out["gene"].astype(str)
    out["log2_fc"] = pd.to_numeric(out["posterior_mean"], errors="coerce") / math.log(2.0)
    out["perturbo_fc_se"] = pd.to_numeric(out["posterior_scale"], errors="coerce") / math.log(2.0)
    return out


def _read_moi(mudata_path: str | Path, override: str | None) -> str | None:
    """The pipeline's multiplicity-of-infection setting, stored by create_mdata."""
    if override:
        return str(override).strip().lower()
    with _open_mudata(mudata_path, backed="r") as mdata:
        raw = mdata[GUIDE_MODALITY].uns.get("moi")
    if raw is None:
        return None
    value = np.asarray(raw).ravel()
    if value.size == 0:
        return None
    text = value[0].decode() if isinstance(value[0], (bytes, np.bytes_)) else str(value[0])
    return text.strip().lower() or None


def _observed_singleton(mudata_path: str | Path) -> tuple[bool, int]:
    """Whether every cell carries at most one perturbation, and how many carry only
    controls. The assignments decide the design; the declared setting can be wrong."""
    with _open_mudata(mudata_path, backed="r") as mdata:
        guide = mdata[GUIDE_MODALITY]
        matrix = guide.layers[ASSIGNMENT_LAYER] if ASSIGNMENT_LAYER in guide.layers else guide.X
        matrix = matrix if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
        per_cell = np.asarray((matrix > 0).sum(axis=1)).ravel()
        control = _guide_control_mask(pd.DataFrame(guide.var)).to_numpy()
        control_only = int(
            (
                (np.asarray((matrix[:, control] > 0).sum(axis=1)).ravel() > 0)
                & (np.asarray((matrix[:, ~control] > 0).sum(axis=1)).ravel() == 0)
            ).sum()
        )
    return bool(per_cell.size and per_cell.max() <= 1), control_only


def _resolve_crt_pool(requested: str, moi: str | None) -> str:
    """Which cells a perturbation is tested against.

    ``from-moi`` keeps the pipeline's own design setting in charge: ``high`` is the
    all-cells pool (every element a marginal association over all cells, the
    behaviour the pipeline has always had), ``low`` the control-anchored pool
    (each perturbation against the pure control cells plus its own). Anything
    else falls back to PerTurbo measuring the design from the data.
    """
    if requested != "from-moi":
        return requested
    if moi == "high":
        return "all-cells"
    if moi == "low":
        return "control-anchored"
    print(f"MOI setting {moi!r} is neither high nor low; letting PerTurbo measure the design (--crt-pool auto).")
    return "auto"


def _record_control_group(artifact_dir: Path, record: dict | None) -> None:
    """Write the control-group provenance beside the PerTurbo artifacts.

    PerTurbo's own ``crt_metadata.json`` already names the pool it used; what it
    cannot know is where that pool came from. Keep that pipeline-specific
    record at the top level; published native JSON stays immutable so its
    recovery manifest remains valid.
    """
    if not record:
        return
    (artifact_dir / "control_group_resolution.json").write_text(
        json.dumps(record, indent=2, default=str)
    )


def _run_perturbo(
    input_path: Path,
    out_dir: Path,
    *,
    map_key: str,
    names_key: str,
    pairs_to_test_path: Path | None,
    gpu_id: str | None,
    phase: str,
    args: argparse.Namespace,
) -> subprocess.Popen:
    cmd = [
        "perturbo",
        "--input",
        str(input_path),
        "--out-dir",
        str(out_dir),
        "--modality-key",
        GENE_MODALITY,
        "--perturbation-modality-key",
        GUIDE_MODALITY,
        "--perturbation-layer",
        ASSIGNMENT_LAYER,
        "--perturbation-element-varm-key",
        map_key,
        "--perturbation-element-names-uns-key",
        names_key,
        "--control-substring",
        CONTROL_SUBSTRING,
        "--library-size-key",
        "total_gene_umis",
        "--size-factor-mode",
        args.size_factor_mode,
        "--likelihood",
        args.likelihood,
        "--prior",
        args.prior,
        "--num-steps-control",
        str(args.num_steps_control),
        "--num-steps-betas",
        str(args.num_steps_betas),
        "--minibatch-size",
        str(args.batch_size),
        "--max-chunk-size",
        str(args.max_chunk_size),
        "--perturbation-chunk-size",
        str(args.perturbation_chunk_size),
        "--device",
        args.device,
        "--step-size",
        str(args.step_size),
        "--no-progress-bar",
    ]
    if args.gene_chunk_size > 0:
        cmd.extend(["--gene-chunk-size", str(args.gene_chunk_size)])
    if args.crt_gene_chunk_size > 0:
        cmd.extend(["--crt-gene-chunk-size", str(args.crt_gene_chunk_size)])
    if args.crt_max_gather_gib > 0:
        cmd.extend(["--crt-max-gather-gib", str(args.crt_max_gather_gib)])
    if pairs_to_test_path is not None:
        # Since PerTurbo 2.0 this does not restrict the fit: it selects the rows of a
        # second table, element_effects_requested_pairs.parquet, with q-values
        # corrected within that set. One run yields the cis-scale and the
        # transcriptome-wide tables together.
        cmd.extend(["--pairs-to-test", str(pairs_to_test_path)])
    if args.crt:
        # The same CRT configuration as the production Gasperini runs: the baseline is
        # polished onto the control null mode by Fisher scoring, and the null-mode guard
        # is advisory rather than fatal. On a full gene panel the guard's percentile is
        # dominated by genes with almost no control counts, which the polish cannot
        # move and the test cannot resolve anyway.
        cmd.extend([
            "--crt", "--crt-mechanism", "propensity", "--crt-tail-families", "saddlepoint",
            "--crt-saddlepoint-only", "--crt-polish-baseline", "--crt-allow-unconverged-baseline",
            "--crt-pool", args.resolved_crt_pool,
        ])
        if args.crt_test_control_elements:
            # The control evaluation scores non-targeting pairs against direct-target
            # pairs, so it needs the controls to carry p-values of their own.
            cmd.append("--crt-test-control-elements")
    # The same covariates SCEPTRE receives, from the same list, so a difference in
    # calls is a difference in method rather than in model specification. The
    # preparation step writes them into the MuData's top-level obs for SCEPTRE and
    # leaves them on the modalities for PerTurbo.
    # Detect on the PerTurbo-side column: for the count terms that is the
    # precomputed `log_` column, not the raw count SCEPTRE reads. `_ensure_covariates`
    # has already guaranteed those exist, including on inputs that skipped the shared
    # preparation step.
    with _open_mudata(input_path, backed="r") as mdata:
        obs_columns = set(mdata[GENE_MODALITY].obs.columns)
    present = [
        c
        for c in CANONICAL_COVARIATES
        if c.perturbo_name is not None and c.perturbo_name in obs_columns
    ]
    continuous, batch = perturbo_covariate_arguments(present)
    # Live now that the guide-UMI term is a conditioned covariate: a covariate
    # constant across the cells PerTurbo fits its baseline on is unidentifiable
    # there. Which cells those are depends on the resolved CRT pool.
    if GUIDE_UMI_COVARIATE in continuous and not _covariate_has_baseline_variance(
        input_path, map_key, names_key, pool=args.resolved_crt_pool
    ):
        scope = (
            "every cell" if args.resolved_crt_pool == "all-cells" else "the control cells"
        )
        print(
            f"Dropping {GUIDE_UMI_COVARIATE}: it has zero variance across {scope}, "
            f"which PerTurbo fits its baseline on under the "
            f"{args.resolved_crt_pool!r} pool."
        )
        continuous = [c for c in continuous if c != GUIDE_UMI_COVARIATE]
    if continuous:
        cmd.extend(["--continuous-covariates", *continuous])
    if batch:
        cmd.extend(["--batch-covariate", batch])
    print(f"Covariates: continuous {continuous or 'none'}; batch {batch or 'none'}.")
    if not args.save_model_params:
        cmd.append("--no-save-model-params")
    out_dir.mkdir(parents=True, exist_ok=True)
    process_env = os.environ.copy()
    if gpu_id is not None:
        process_env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    cache_root = Path(args.jax_cache_dir)
    if not cache_root.is_absolute():
        cache_root = Path.cwd() / cache_root
    cache_dir = cache_root / phase
    cache_dir.mkdir(parents=True, exist_ok=True)
    process_env["JAX_COMPILATION_CACHE_DIR"] = str(cache_dir)
    process_env["JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS"] = "0"
    process_env["JAX_PERSISTENT_CACHE_MIN_ENTRY_SIZE_BYTES"] = "-1"
    process_env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    phase_tmp = out_dir / "tmp"
    phase_tmp.mkdir(parents=True, exist_ok=True)
    process_env["TMPDIR"] = str(phase_tmp)
    print(
        f"Running PerTurbo v2 {phase} fit on CUDA_VISIBLE_DEVICES={process_env.get('CUDA_VISIBLE_DEVICES', '<inherited>')}: "
        + " ".join(cmd)
    )
    return subprocess.Popen(cmd, env=process_env)


def _wait_for_fits(processes: dict[str, subprocess.Popen]) -> None:
    pending = dict(processes)
    while pending:
        for name, process in list(pending.items()):
            return_code = process.poll()
            if return_code is None:
                continue
            pending.pop(name)
            if return_code != 0:
                for peer in pending.values():
                    peer.terminate()
                for peer in pending.values():
                    peer.wait()
                raise subprocess.CalledProcessError(return_code, process.args)
            print(f"PerTurbo v2 {name} fit completed successfully.")
        if pending:
            time.sleep(1)


def _requested_pairs_table(fit_dir: Path, effects: pd.DataFrame, pairs_path: Path | None) -> pd.DataFrame:
    """PerTurbo's own requested-pairs table, whose q-values are corrected within the
    requested set; falls back to selecting the rows here when an older PerTurbo did
    not write it (convert_* recomputes the q-values from p either way)."""
    written = fit_dir / "element_effects_requested_pairs.parquet"
    if written.exists():
        return pd.read_parquet(written)
    if pairs_path is None:
        return effects
    pairs = pd.read_parquet(pairs_path)[["element", "gene"]].astype(str).drop_duplicates()
    keyed = effects.assign(element=effects["element"].astype(str), gene=effects["gene"].astype(str))
    return keyed.merge(pairs, on=["element", "gene"], how="inner")


def _pvalue_projection(raw_path: Path, output_path: Path, *, crt: bool, batch_rows: int) -> Path:
    """Project only the selected p-value, in bounded Arrow batches."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    parquet = pq.ParquetFile(raw_path)
    names = set(parquet.schema_arrow.names)
    candidates = ["crt_saddlepoint_p_value", "crt_p_value", "empirical_p_value"]
    selected = next((name for name in candidates if name in names), None)
    columns = [name for name in (selected, "posterior_prob") if name is not None and name in names]
    temporary = output_path.with_name(output_path.name + ".partial")
    writer = None
    try:
        for batch in parquet.iter_batches(batch_size=batch_rows, columns=columns):
            frame = batch.to_pandas()
            p = pd.to_numeric(frame[selected], errors="coerce") if selected else pd.Series(np.nan, index=frame.index, dtype=float)
            if not crt:
                posterior = pd.to_numeric(frame["posterior_prob"], errors="coerce")
                p = p.where(p.notna(), posterior)
            # Match the target branch's legacy _bh_adjust exactly: every non-NaN
            # value, including +/-inf, participates after clipping. The unmodified
            # p-value is still written to the result table in the second pass.
            projected = pd.DataFrame({"p_value": p.clip(0.0, 1.0).astype("float64", copy=False)})
            table = pa.Table.from_pandas(projected, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(temporary, table.schema, compression="zstd")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        pd.DataFrame({"p_value": pd.Series(dtype="float64")}).to_parquet(temporary, index=False)
    os.replace(temporary, output_path)
    return output_path


def _bounded_raw_columns(parquet) -> list[str]:
    names = set(parquet.schema_arrow.names)
    required = ["element", "gene", "posterior_mean", "posterior_scale", "posterior_prob"]
    missing = sorted(set(required) - names)
    if missing:
        raise ValueError("PerTurbo result is missing columns: " + ", ".join(missing))
    p_column = next(
        (name for name in ("crt_saddlepoint_p_value", "crt_p_value", "empirical_p_value") if name in names),
        None,
    )
    diagnostics = [name for name in CRT_DIAGNOSTIC_COLUMNS if name in names]
    return required + ([p_column] if p_column else []) + diagnostics


def _ensure_requested_pairs(raw_path: Path, requested_path: Path, pairs_path: Path | None, *, batch_rows: int) -> Path:
    """Build the older-PerTurbo fallback with a bounded join to the small pair set."""
    if requested_path.exists():
        return requested_path
    if pairs_path is None:
        raise FileNotFoundError("Requested-pairs output is absent and no pairs table was supplied")
    import pyarrow as pa
    import pyarrow.parquet as pq

    pairs = pd.read_parquet(pairs_path, columns=["element", "gene"]).astype(str).drop_duplicates()
    parquet = pq.ParquetFile(raw_path)
    columns = _bounded_raw_columns(parquet)
    temporary = requested_path.with_name(requested_path.name + ".partial")
    writer = None
    try:
        for batch in parquet.iter_batches(batch_size=batch_rows, columns=columns):
            frame = batch.to_pandas()
            frame["element"] = frame["element"].astype(str)
            frame["gene"] = frame["gene"].astype(str)
            selected = frame.merge(pairs, on=["element", "gene"], how="inner")
            table = pa.Table.from_pandas(selected, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(temporary, table.schema, compression="zstd")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        pd.DataFrame(columns=columns).to_parquet(temporary, index=False)
    os.replace(temporary, requested_path)
    return requested_path


def _write_bounded_effects(
    raw_path: Path,
    output_path: str | Path,
    *,
    inference_type: str,
    prepared_path: Path | None,
    guide_name_map: dict[str, str],
    crt: bool,
    scratch_dir: Path,
    max_bh_working_bytes: int,
    compact_floats: bool,
    batch_rows: int,
    target_lookup: pd.DataFrame | None = None,
) -> None:
    """Convert one complete family without materializing its wide raw table."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    started = time.perf_counter()
    output_path = Path(output_path)
    print(f"[perturbo conversion] {inference_type}: {raw_path} -> {output_path}; batch_rows={batch_rows}", flush=True)
    p_path = _pvalue_projection(raw_path, scratch_dir / f"{inference_type}.p.parquet", crt=crt, batch_rows=batch_rows)
    q_path = write_bh_sidecar(
        p_path,
        scratch_dir / f"{inference_type}.q.parquet",
        max_working_bytes=max_bh_working_bytes,
    )
    q_values = pq.read_table(q_path, columns=["perturbo_q_value"]).column(0).to_numpy()
    if inference_type == "element":
        if target_lookup is None:
            if prepared_path is None:
                raise ValueError("Element conversion needs target_lookup or prepared_path")
            with _open_mudata(prepared_path, backed="r") as mdata:
                target_lookup = get_target_lookup(pd.DataFrame(mdata[GUIDE_MODALITY].var).copy())
        else:
            target_lookup = target_lookup.copy()
        if target_lookup["intended_target_key"].duplicated().any():
            raise ValueError(
                "Bounded conversion requires one metadata row per intended_target_key; "
                "duplicates would change the BH family during the legacy merge."
            )
        output_columns = ["gene_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end", "log2_fc", "perturbo_fc_se", "p_value", "perturbo_posterior_prob"]
    else:
        output_columns = ["gene_id", "guide_id", "log2_fc", "perturbo_fc_se", "p_value", "perturbo_posterior_prob"]

    raw_schema = pq.ParquetFile(raw_path).schema_arrow
    present_diagnostics = [name for name in CRT_DIAGNOSTIC_COLUMNS if name in raw_schema.names]
    output_columns.extend(f"perturbo_{name}" for name in present_diagnostics)
    output_columns.append("perturbo_q_value")

    def _result_float_type(raw_name: str, *, compact: bool = False):
        raw_type = raw_schema.field(raw_name).type
        if compact:
            return pa.float32()
        return pa.float32() if pa.types.is_float32(raw_type) else pa.float64()

    dictionary_string = pa.dictionary(pa.int32(), pa.string())
    fields = [pa.field("gene_id", dictionary_string)]
    if inference_type == "element":
        fields.extend([
            pa.field("intended_target_name", dictionary_string),
            pa.field("intended_target_chr", dictionary_string),
        ])
        lookup_schema = pa.Schema.from_pandas(
            target_lookup[["intended_target_start", "intended_target_end"]],
            preserve_index=False,
        )
        for coordinate in ("intended_target_start", "intended_target_end"):
            coordinate_type = lookup_schema.field(coordinate).type
            if pa.types.is_null(coordinate_type):
                coordinate_type = pa.float64()
            fields.append(pa.field(coordinate, coordinate_type))
    else:
        fields.append(pa.field("guide_id", dictionary_string))
    fields.extend([
        pa.field("log2_fc", _result_float_type("posterior_mean", compact=compact_floats)),
        pa.field("perturbo_fc_se", _result_float_type("posterior_scale", compact=compact_floats)),
        pa.field("p_value", pa.float64()),
        pa.field("perturbo_posterior_prob", _result_float_type("posterior_prob")),
    ])
    fields.extend(
        pa.field(f"perturbo_{name}", raw_schema.field(name).type)
        for name in present_diagnostics
    )
    fields.append(pa.field("perturbo_q_value", pa.float64()))
    output_schema = pa.schema(fields)

    def _table_with_output_schema(frame: pd.DataFrame):
        arrays = []
        for field in output_schema:
            if pa.types.is_dictionary(field.type):
                values = pa.array(frame[field.name], type=field.type.value_type, from_pandas=True)
                arrays.append(values.dictionary_encode())
            else:
                arrays.append(pa.array(frame[field.name], type=field.type, from_pandas=True))
        return pa.Table.from_arrays(arrays, schema=output_schema)

    temporary = output_path.with_name(output_path.name + ".partial")
    parquet_writer = None
    text_handle = None
    offset = 0
    try:
        if str(output_path).lower().endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            text_handle = gzip.open(temporary, "wt") if str(output_path).lower().endswith(".gz") else temporary.open("w")
        raw_parquet = pq.ParquetFile(raw_path)
        raw_columns = _bounded_raw_columns(raw_parquet)
        for batch in raw_parquet.iter_batches(batch_size=batch_rows, columns=raw_columns):
            frame = _convert_common_effect_columns(batch.to_pandas(), crt=crt)
            n = len(frame)
            frame["perturbo_q_value"] = q_values[offset : offset + n]
            offset += n
            if inference_type == "element":
                frame = frame.rename(columns={"element": "intended_target_key"}).merge(
                    target_lookup, on="intended_target_key", how="left"
                )
                if frame["intended_target_name"].isna().any():
                    raise ValueError("Unable to map one or more PerTurbo element names")
            else:
                frame["guide_id"] = frame["element"].map(guide_name_map).fillna(frame["element"]).astype(str)
            frame = frame[output_columns]
            if compact_floats:
                frame = compact_result_floats(frame)
            if text_handle is not None:
                frame.to_csv(text_handle, sep="\t", index=False, header=(offset == n))
            else:
                table = _table_with_output_schema(frame)
                if parquet_writer is None:
                    parquet_writer = pq.ParquetWriter(temporary, output_schema, compression="zstd")
                parquet_writer.write_table(table)
    finally:
        if parquet_writer is not None:
            parquet_writer.close()
        if text_handle is not None:
            text_handle.close()
    if offset != len(q_values):
        raise ValueError("Raw result and BH sidecar row counts differ")
    if offset == 0:
        if text_handle is None:
            pq.write_table(_table_with_output_schema(pd.DataFrame(columns=output_columns)), temporary)
        else:
            empty = pd.DataFrame(columns=output_columns)
            with (gzip.open(temporary, "wt") if str(output_path).lower().endswith(".gz") else temporary.open("w")) as handle:
                empty.to_csv(handle, sep="\t", index=False)
    os.replace(temporary, output_path)
    import resource
    import sys

    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mib = peak_rss / (1024**2 if sys.platform == "darwin" else 1024)
    print(
        f"[perturbo conversion] {inference_type}: {offset:,} rows in "
        f"{time.perf_counter() - started:.2f}s; adapter_lifetime_peak_rss_mib={peak_mib:.1f}",
        flush=True,
    )


def _publish_raw_fit(source: Path, artifact_root: Path, name: str) -> Path:
    """Durably publish one successful fit before another fit can start."""
    target = artifact_root / name
    partial = artifact_root / f".{name}.partial"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite durable PerTurbo fit: {target}")
    if source.resolve() != partial.resolve():
        if partial.exists():
            shutil.rmtree(partial)
        shutil.copytree(source, partial)
    manifest = {"raw_parquets": [parquet_manifest(partial / "element_effects.parquet")]}
    requested = partial / "element_effects_requested_pairs.parquet"
    if requested.exists():
        manifest["raw_parquets"].append(parquet_manifest(requested))
    manifest["source_metadata"] = [
        {"path": str(path.relative_to(partial)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for path in sorted(partial.rglob("*.json"))
    ]
    write_manifest(partial / "conversion_manifest.json", manifest)
    os.replace(partial, target)
    return target


def _verify_durable_fit(directory: Path) -> None:
    manifest = json.loads((directory / "conversion_manifest.json").read_text())
    verify_manifest(directory, manifest)
    recorded_parquets = {entry["path"] for entry in manifest["raw_parquets"]}
    requested = directory / "element_effects_requested_pairs.parquet"
    if requested.exists() and requested.name not in recorded_parquets:
        raise ValueError(f"Unverified requested-pairs artifact in durable fit: {requested}")
    for expected in manifest.get("source_metadata", []):
        path = directory / expected["path"]
        observed = {"path": expected["path"], "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        if observed != expected:
            raise ValueError(f"PerTurbo source metadata provenance failed: {path}")


def _atomic_copy(source: Path, destination: str | Path) -> None:
    destination = Path(destination)
    if source.resolve() == destination.resolve():
        return
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copy2(source, temporary)
    os.replace(temporary, destination)


def _result_encoding(path: str | Path) -> str:
    value = str(path).lower()
    if value.endswith((".parquet", ".pq")):
        return "parquet"
    return "tsv.gz" if value.endswith(".gz") else "tsv"


def _wait_publish_fits(processes, fit_dirs, artifact_dir):
    """Publish every successful fit immediately, even when a peer later fails."""
    pending = dict(processes)
    published = {}
    try:
        while pending:
            completed = []
            failures = []
            for name, process in list(pending.items()):
                return_code = process.poll()
                if return_code is None:
                    continue
                completed.append(name)
                if return_code == 0:
                    published[name] = _publish_raw_fit(fit_dirs[name], artifact_dir, name)
                    print(f"PerTurbo v2 {name} fit completed and was durably published.")
                else:
                    failures.append((name, process, return_code))
            for name in completed:
                pending.pop(name)
            if failures:
                # One last poll closes the race in which a peer succeeded while the
                # failing process was being handled; preserve that completed work.
                for peer_name, peer in list(pending.items()):
                    if peer.poll() == 0:
                        published[peer_name] = _publish_raw_fit(
                            fit_dirs[peer_name], artifact_dir, peer_name
                        )
                        pending.pop(peer_name)
                for peer in pending.values():
                    peer.terminate()
                for peer in pending.values():
                    peer.wait()
                _, process, return_code = failures[0]
                raise subprocess.CalledProcessError(return_code, process.args)
            if pending:
                time.sleep(1)
    except BaseException:
        for peer in pending.values():
            if peer.poll() is None:
                peer.terminate()
        for peer in pending.values():
            peer.wait()
        raise
    return published


def _fit_provenance(args, guide_name_map, element_pairs_path, guide_pairs_path) -> dict:
    source = Path(args.input).resolve()
    stat = source.stat()
    return {
        "input": {"path": str(source), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns},
        "guide_name_map": sorted([str(k), str(v)] for k, v in guide_name_map.items()),
        "element_pairs_sha256": sha256_file(element_pairs_path) if element_pairs_path else None,
        "guide_pairs_sha256": sha256_file(guide_pairs_path) if guide_pairs_path else None,
        "fit_settings": {
            name: getattr(args, name)
            for name in (
                "crt", "resolved_crt_pool", "size_factor_mode", "likelihood", "prior",
                "num_steps_control", "num_steps_betas", "step_size", "max_chunk_size",
                "perturbation_chunk_size", "gene_chunk_size", "crt_gene_chunk_size",
                "crt_max_gather_gib", "batch_size",
            )
        },
    }


def _write_conversion_metadata(artifact_dir, prepared, guide_name_map, element_pairs_path, guide_pairs_path):
    """Persist the compact mappings needed to convert without reopening MuData."""
    metadata_dir = artifact_dir / "conversion_metadata"
    metadata_dir.mkdir(parents=True, exist_ok=False)
    with _open_mudata(prepared, backed="r") as mdata:
        target_lookup = get_target_lookup(pd.DataFrame(mdata[GUIDE_MODALITY].var).copy())
    target_path = metadata_dir / "target_lookup.parquet"
    target_lookup.to_parquet(target_path, index=False)
    guide_path = metadata_dir / "guide_name_map.json"
    guide_path.write_text(json.dumps(guide_name_map, sort_keys=True) + "\n")
    paths = [target_path, guide_path]
    for name, source in (
        ("element_pairs_to_test.parquet", element_pairs_path),
        ("guide_pairs_to_test.parquet", guide_pairs_path),
    ):
        if source is not None:
            destination = metadata_dir / name
            shutil.copy2(source, destination)
            paths.append(destination)
    write_manifest(
        metadata_dir / "manifest.json",
        {"files": [{"path": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in paths]},
    )
    return metadata_dir


def _load_conversion_metadata(artifact_dir):
    metadata_dir = artifact_dir / "conversion_metadata"
    manifest_path = metadata_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            "Durable artifacts predate compact conversion metadata and cannot be "
            "reinterpreted automatically; rerun the fit or adopt them manually."
        )
    manifest = json.loads(manifest_path.read_text())
    for expected in manifest["files"]:
        path = metadata_dir / expected["path"]
        observed = {"path": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        if observed != expected:
            raise ValueError(f"Conversion metadata provenance failed: {path}")
    return {
        "guide_name_map": json.loads((metadata_dir / "guide_name_map.json").read_text()),
        "target_lookup": pd.read_parquet(metadata_dir / "target_lookup.parquet"),
        "element_pairs": (metadata_dir / "element_pairs_to_test.parquet") if (metadata_dir / "element_pairs_to_test.parquet").exists() else None,
        "guide_pairs": (metadata_dir / "guide_pairs_to_test.parquet") if (metadata_dir / "guide_pairs_to_test.parquet").exists() else None,
    }


def _run_conversion_only(args) -> None:
    """Rebuild tables from durable compact metadata without loading the input MuData."""
    if args.output_mudata:
        raise ValueError("Conversion-only recovery cannot write --output-mudata without loading the full input")
    artifact_dir = Path(args.conversion_only_artifact_dir)
    recorded = json.loads((artifact_dir / "fit_provenance.json").read_text())
    source = Path(args.input).resolve()
    stat = source.stat()
    observed_input = {"path": str(source), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if observed_input != recorded["input"]:
        raise ValueError("Conversion-only input provenance does not match the producing fit")
    for name, value in recorded["fit_settings"].items():
        if name == "resolved_crt_pool":
            continue
        if getattr(args, name) != value:
            raise ValueError(f"Conversion-only fit setting changed: {name}")
    args.resolved_crt_pool = recorded["fit_settings"]["resolved_crt_pool"]
    metadata = _load_conversion_metadata(artifact_dir)
    durable_element_dir = artifact_dir / "element"
    durable_guide_dir = artifact_dir / "guide"
    for durable in (durable_element_dir, durable_guide_dir):
        _verify_durable_fit(durable)
    wants_local = (not args.test_all_pairs) or bool(args.local_per_element_output) or bool(args.local_per_guide_output)
    with tempfile.TemporaryDirectory(prefix="perturbo_v2_conversion_") as tmp:
        scratch = Path(tmp)
        if args.test_all_pairs:
            _write_bounded_effects(durable_element_dir / "element_effects.parquet", Path(args.per_element_output), inference_type="element", prepared_path=None, target_lookup=metadata["target_lookup"], guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            _write_bounded_effects(durable_guide_dir / "element_effects.parquet", Path(args.per_guide_output), inference_type="guide", prepared_path=None, guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
        if wants_local:
            native_element_requested = durable_element_dir / "element_effects_requested_pairs.parquet"
            native_guide_requested = durable_guide_dir / "element_effects_requested_pairs.parquet"
            element_requested = native_element_requested if native_element_requested.exists() else _ensure_requested_pairs(durable_element_dir / "element_effects.parquet", scratch / "element_requested.parquet", metadata["element_pairs"], batch_rows=args.result_batch_rows)
            guide_requested = native_guide_requested if native_guide_requested.exists() else _ensure_requested_pairs(durable_guide_dir / "element_effects.parquet", scratch / "guide_requested.parquet", metadata["guide_pairs"], batch_rows=args.result_batch_rows)
            element_output = Path(args.per_element_output) if not args.test_all_pairs else Path(args.local_per_element_output or (scratch / "local_element.parquet"))
            guide_output = Path(args.per_guide_output) if not args.test_all_pairs else Path(args.local_per_guide_output or (scratch / "local_guide.parquet"))
            _write_bounded_effects(element_requested, element_output, inference_type="element", prepared_path=None, target_lookup=metadata["target_lookup"], guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            _write_bounded_effects(guide_requested, guide_output, inference_type="guide", prepared_path=None, guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            if not args.test_all_pairs:
                if args.local_per_element_output:
                    if _result_encoding(element_output) == _result_encoding(args.local_per_element_output):
                        _atomic_copy(element_output, args.local_per_element_output)
                    else:
                        _write_bounded_effects(element_requested, args.local_per_element_output, inference_type="element", prepared_path=None, target_lookup=metadata["target_lookup"], guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
                if args.local_per_guide_output:
                    if _result_encoding(guide_output) == _result_encoding(args.local_per_guide_output):
                        _atomic_copy(guide_output, args.local_per_guide_output)
                    else:
                        _write_bounded_effects(guide_requested, args.local_per_guide_output, inference_type="guide", prepared_path=None, guide_name_map=metadata["guide_name_map"], crt=args.crt, scratch_dir=scratch, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)

    write_manifest(artifact_dir / "last_conversion.json", {
        "conversion_only": True,
        "compact_result_floats": args.compact_result_floats,
        "max_bh_working_bytes": args.max_bh_working_bytes,
        "result_batch_rows": args.result_batch_rows,
    })


def run_pipeline_adapter(args: argparse.Namespace) -> None:
    if args.result_batch_rows <= 0:
        raise ValueError("--result-batch-rows must be positive")
    if args.max_bh_working_bytes <= 0:
        raise ValueError("--max-bh-working-bytes must be positive")
    if args.conversion_only_artifact_dir:
        _run_conversion_only(args)
        return
    with tempfile.TemporaryDirectory(prefix="perturbo_v2_pipeline_") as tmp:
        tmp_dir = Path(tmp)
        prepared = tmp_dir / "prepared_mudata.h5mu"
        guide_name_map = prepare_mudata_for_perturbo_v2(
            args.input,
            prepared,
            test_all_pairs=args.test_all_pairs,
        )

        declared_moi = _read_moi(args.input, args.moi)
        args.resolved_crt_pool = _resolve_crt_pool(args.crt_pool, declared_moi)
        pool_reason = "requested outright" if args.crt_pool != "from-moi" else f"mapped from declared MOI {declared_moi!r}"
        if args.crt_pool == "from-moi":
            # The assignments outrank the declared setting. A screen whose cells each
            # carry one perturbation is a low-MOI screen whatever the samplesheet says,
            # and SCEPTRE's driver reaches the same conclusion from the same matrix, so
            # both methods contrast against the control cells rather than one against
            # the complement.
            singleton, control_only = _observed_singleton(args.input)
            if singleton and control_only > 0 and args.resolved_crt_pool != "control-anchored":
                print(
                    f"Every cell carries at most one perturbation and {control_only:,} carry only "
                    f"controls: using the control-anchored pool rather than "
                    f"'{args.resolved_crt_pool}' from the declared MOI."
                )
                args.resolved_crt_pool = "control-anchored"
                pool_reason = (
                    f"every cell carries at most one perturbation and {control_only:,} carry only controls"
                )
        # The same facts the pipeline logged, written beside the results so a run can
        # be read back without its log: which cells the test compared against, the
        # shared setting it came from, and whether that setting was taken from the
        # MOI, stated outright, or overridden for PerTurbo alone.
        args.control_group_record = {
            "method": "perturbo",
            "crt_pool": args.resolved_crt_pool,
            "crt_pool_requested": args.crt_pool,
            "declared_moi": normalize_moi(declared_moi),
            SETTING_PARAM: args.control_group_setting,
            "provenance": args.control_group_provenance,
            "reason": pool_reason,
        }
        if args.crt:
            print(f"PerTurbo CRT pool: {args.resolved_crt_pool} (requested {args.crt_pool}; {pool_reason}).")
        # The requested pairs (a cis window, usually) come from the MuData that
        # carries uns['pairs_to_test']; by default that is the input itself. The fit
        # is never restricted to them - they select the rows of the local tables.
        element_pairs_path: Path | None = None
        guide_pairs_path: Path | None = None
        wants_local = (not args.test_all_pairs) or bool(args.local_per_element_output) or bool(args.local_per_guide_output)
        if wants_local:
            pairs_source = args.pairs_mudata or args.input
            with _open_mudata(pairs_source, backed="r") as pairs_mdata:
                pairs_frame = pairs_mdata.uns.get("pairs_to_test")
                if pairs_frame is None:
                    raise KeyError(f"pairs_to_test not found in {pairs_source}; the local tables need it.")
                pairs_frame = pairs_frame.copy() if isinstance(pairs_frame, pd.DataFrame) else pd.DataFrame(pairs_frame)
            with _open_mudata(prepared) as prepared_mdata:
                prepared_mdata.uns["pairs_to_test"] = pairs_frame
                element_pairs_path = tmp_dir / "element_pairs_to_test.parquet"
                guide_pairs_path = tmp_dir / "guide_pairs_to_test.parquet"
                _build_native_pairs_to_test(
                    prepared_mdata,
                    pair_element_column="intended_target_key",
                    element_names=[str(x) for x in prepared_mdata[GUIDE_MODALITY].uns[ELEMENT_NAMES_KEY]],
                ).to_parquet(element_pairs_path, index=False)
                _build_native_pairs_to_test(
                    prepared_mdata,
                    pair_element_column="guide_id",
                    element_names=[str(x) for x in prepared_mdata[GUIDE_MODALITY].uns[GUIDE_NAMES_KEY]],
                    guide_name_map=guide_name_map,
                ).to_parquet(guide_pairs_path, index=False)

        artifact_dir = Path(args.conversion_only_artifact_dir or args.v2_artifact_dir or (Path(args.per_element_output).parent / "perturbo_v2_outputs"))
        fit_provenance = _fit_provenance(args, guide_name_map, element_pairs_path, guide_pairs_path)
        if args.conversion_only_artifact_dir:
            recorded_provenance = json.loads((artifact_dir / "fit_provenance.json").read_text())
            if recorded_provenance != fit_provenance:
                raise ValueError("Conversion-only fit provenance does not match input, mappings, pairs, or fit settings")
            durable_element_dir = artifact_dir / "element"
            durable_guide_dir = artifact_dir / "guide"
            for durable in (durable_element_dir, durable_guide_dir):
                _verify_durable_fit(durable)
            print(f"Verified durable PerTurbo artifacts for conversion-only recovery: {artifact_dir}")
        else:
            artifact_dir.mkdir(parents=True, exist_ok=True)
            reserved = [artifact_dir / name for name in ("element", "guide", ".element.partial", ".guide.partial", "fit_provenance.json")]
            collisions = [path for path in reserved if path.exists()]
            if collisions:
                raise FileExistsError(
                    "Refusing to overwrite an existing PerTurbo fit; use "
                    f"--conversion-only-artifact-dir for recovery: {collisions[0]}"
                )
            write_manifest(artifact_dir / "fit_provenance.json", fit_provenance)
            _write_conversion_metadata(
                artifact_dir, prepared, guide_name_map,
                element_pairs_path, guide_pairs_path,
            )
            element_dir = artifact_dir / ".element.partial"
            guide_dir = artifact_dir / ".guide.partial"
            element_process = _run_perturbo(
                prepared,
                element_dir,
                map_key=ELEMENT_MAP_KEY,
                names_key=ELEMENT_NAMES_KEY,
                pairs_to_test_path=element_pairs_path,
                gpu_id=args.element_gpu if args.parallel_fits else None,
                phase="element",
                args=args,
            )
            if args.parallel_fits:
                guide_process = _run_perturbo(
                    prepared, guide_dir, map_key=GUIDE_MAP_KEY, names_key=GUIDE_NAMES_KEY,
                    pairs_to_test_path=guide_pairs_path, gpu_id=args.guide_gpu,
                    phase="guide", args=args,
                )
                published = _wait_publish_fits(
                    {"element": element_process, "guide": guide_process},
                    {"element": element_dir, "guide": guide_dir}, artifact_dir,
                )
                durable_element_dir = published["element"]
                durable_guide_dir = published["guide"]
            else:
                durable_element_dir = _wait_publish_fits(
                    {"element": element_process}, {"element": element_dir}, artifact_dir
                )["element"]
                guide_process = _run_perturbo(
                    prepared, guide_dir, map_key=GUIDE_MAP_KEY, names_key=GUIDE_NAMES_KEY,
                    pairs_to_test_path=guide_pairs_path, gpu_id=None, phase="guide", args=args,
                )
                durable_guide_dir = _wait_publish_fits(
                    {"guide": guide_process}, {"guide": guide_dir}, artifact_dir
                )["guide"]

        global_element_path = global_guide_path = None
        if args.test_all_pairs or args.output_mudata:
            global_element_path = tmp_dir / "global_element.parquet" if args.output_mudata else Path(args.per_element_output)
            global_guide_path = tmp_dir / "global_guide.parquet" if args.output_mudata else Path(args.per_guide_output)
            _write_bounded_effects(durable_element_dir / "element_effects.parquet", global_element_path, inference_type="element", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            _write_bounded_effects(durable_guide_dir / "element_effects.parquet", global_guide_path, inference_type="guide", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            if args.test_all_pairs and args.output_mudata:
                for raw, canonical, destination, kind in (
                    (durable_element_dir / "element_effects.parquet", global_element_path, args.per_element_output, "element"),
                    (durable_guide_dir / "element_effects.parquet", global_guide_path, args.per_guide_output, "guide"),
                ):
                    if _result_encoding(canonical) == _result_encoding(destination):
                        _atomic_copy(canonical, destination)
                    else:
                        _write_bounded_effects(raw, destination, inference_type=kind, prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)

        local_element_path = local_guide_path = None
        if wants_local:
            native_element_requested = durable_element_dir / "element_effects_requested_pairs.parquet"
            native_guide_requested = durable_guide_dir / "element_effects_requested_pairs.parquet"
            element_requested = native_element_requested if native_element_requested.exists() else _ensure_requested_pairs(
                durable_element_dir / "element_effects.parquet", tmp_dir / "element_requested.parquet",
                element_pairs_path, batch_rows=args.result_batch_rows,
            )
            guide_requested = native_guide_requested if native_guide_requested.exists() else _ensure_requested_pairs(
                durable_guide_dir / "element_effects.parquet", tmp_dir / "guide_requested.parquet",
                guide_pairs_path, batch_rows=args.result_batch_rows,
            )
            local_element_path = tmp_dir / "local_element.parquet" if args.output_mudata else (Path(args.per_element_output) if not args.test_all_pairs else Path(args.local_per_element_output or (tmp_dir / "local_element.parquet")))
            local_guide_path = tmp_dir / "local_guide.parquet" if args.output_mudata else (Path(args.per_guide_output) if not args.test_all_pairs else Path(args.local_per_guide_output or (tmp_dir / "local_guide.parquet")))
            _write_bounded_effects(element_requested, local_element_path, inference_type="element", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            _write_bounded_effects(guide_requested, local_guide_path, inference_type="guide", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            if not args.test_all_pairs and args.output_mudata:
                for raw, canonical, destination, kind in (
                    (element_requested, local_element_path, args.per_element_output, "element"),
                    (guide_requested, local_guide_path, args.per_guide_output, "guide"),
                ):
                    if _result_encoding(canonical) == _result_encoding(destination):
                        _atomic_copy(canonical, destination)
                    else:
                        _write_bounded_effects(raw, destination, inference_type=kind, prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
            if not args.test_all_pairs or args.output_mudata:
                if args.local_per_element_output:
                    if _result_encoding(local_element_path) == _result_encoding(args.local_per_element_output):
                        _atomic_copy(local_element_path, args.local_per_element_output)
                    else:
                        _write_bounded_effects(element_requested, args.local_per_element_output, inference_type="element", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
                if args.local_per_guide_output:
                    if _result_encoding(local_guide_path) == _result_encoding(args.local_per_guide_output):
                        _atomic_copy(local_guide_path, args.local_per_guide_output)
                    else:
                        _write_bounded_effects(guide_requested, args.local_per_guide_output, inference_type="guide", prepared_path=prepared, guide_name_map=guide_name_map, crt=args.crt, scratch_dir=tmp_dir, max_bh_working_bytes=args.max_bh_working_bytes, compact_floats=args.compact_result_floats, batch_rows=args.result_batch_rows)
        write_manifest(
            artifact_dir / "last_conversion.json",
            {
                "compact_result_floats": args.compact_result_floats,
                "max_bh_working_bytes": args.max_bh_working_bytes,
                "result_batch_rows": args.result_batch_rows,
            },
        )
        # --per-*-output keep their historical meaning: global with --test-all-pairs,
        # local without. The --local-per-*-output files add the local tables beside
        # the global ones so a single run serves both.
        if not args.output_mudata:
            print(
                "No --output-mudata given: the result tables are the published output and no MuData "
                "is copied. The pipeline's merge step assembles one from these tables."
            )
        if args.output_mudata:
            print("--output-mudata explicitly requests materializing converted tables in host memory.")
            global_element_df = pd.read_parquet(global_element_path)
            global_guide_df = pd.read_parquet(global_guide_path)
            local_element_df = local_guide_df = None
            if local_element_path is not None:
                local_element_df = pd.read_parquet(local_element_path)
                local_guide_df = pd.read_parquet(local_guide_path)
            # The analysis-qualified keys carry the tables; the generic keys that
            # standalone consumers read are hard links to whichever pair is primary.
            # They are the same frames, and at screen scale a duplicate is gigabytes.
            updates = {
                "global_analysis_per_element_results": make_h5mu_safe_dataframe(global_element_df),
                "global_analysis_per_guide_results": make_h5mu_safe_dataframe(global_guide_df),
            }
            if local_element_df is not None:
                updates["local_analysis_per_element_results"] = make_h5mu_safe_dataframe(local_element_df)
                updates["local_analysis_per_guide_results"] = make_h5mu_safe_dataframe(local_guide_df)
            prefix = "global_analysis" if args.test_all_pairs else "local_analysis"
            aliases = {
                "per_element_results": f"{prefix}_per_element_results",
                "per_guide_results": f"{prefix}_per_guide_results",
            }
            write_uns_patch(args.input, args.output_mudata, updates=updates, aliases=aliases)

        if args.v2_artifact_dir:
            for name, source in {
                "element_pairs_to_test.parquet": element_pairs_path,
                "guide_pairs_to_test.parquet": guide_pairs_path,
            }.items():
                if source is not None:
                    shutil.copy2(source, artifact_dir / name)
            _record_control_group(artifact_dir, getattr(args, "control_group_record", None))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run PerTurbo v2 and emit CRISPR_Pipeline-compatible result files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input CRISPR_Pipeline MuData file")
    parser.add_argument("--per-element-output", required=True, help="Output per-element table (.tsv.gz or .parquet, by extension)")
    parser.add_argument("--per-guide-output", required=True, help="Output per-guide table (.tsv.gz or .parquet, by extension)")
    parser.add_argument(
        "--output-mudata",
        help=(
            "Optional MuData with the result tables in .uns. Writing it byte-copies the input and "
            "re-serialises the tables, which at screen scale is tens of gigabytes; the pipeline "
            "assembles a MuData again downstream, so skip this unless an intermediate consumer "
            "needs one."
        ),
    )
    parser.add_argument("--v2-artifact-dir", default=None, help="Optional directory for raw PerTurbo v2 artifacts")
    parser.add_argument(
        "--conversion-only-artifact-dir",
        default=None,
        help="Skip fitting and rebuild outputs from hash-verified durable raw artifacts",
    )
    parser.add_argument(
        "--max-bh-working-bytes",
        type=int,
        default=24 << 30,
        help="Planning budget for numeric-only family-wide BH working arrays",
    )
    parser.add_argument(
        "--result-batch-rows",
        type=int,
        default=250_000,
        help="Rows per Arrow/Pandas conversion batch",
    )
    parser.add_argument(
        "--compact-result-floats",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Store allowlisted effect/SE fields as float32 after inference and BH",
    )
    parser.add_argument(
        "--test-all-pairs",
        action="store_true",
        help=(
            "Write the transcriptome-wide tables to --per-*-output. The fit always covers every pair; "
            "without this flag --per-*-output hold the requested pairs (mdata.uns['pairs_to_test'])."
        ),
    )
    parser.add_argument("--pairs-mudata", default=None, help="MuData carrying uns['pairs_to_test'] for the local tables (default: --input)")
    parser.add_argument("--local-per-element-output", default=None, help="Also write the requested-pairs per-element table here")
    parser.add_argument("--local-per-guide-output", default=None, help="Also write the requested-pairs per-guide table here")
    parser.add_argument("--crt", action=argparse.BooleanOptionalAction, default=True, help="Run the conditional randomization test (saddlepoint approximation, no resampling)")
    parser.add_argument(
        "--crt-pool",
        default="from-moi",
        choices=["from-moi", "auto", "all-cells", "control-anchored"],
        help="Cells a perturbation is tested against; from-moi maps the pipeline's MOI setting (high -> all-cells, low -> control-anchored)",
    )
    parser.add_argument("--moi", default=None, help="Override the MOI setting read from guide.uns['moi']")
    parser.add_argument(
        "--control-group-setting",
        default=None,
        help=(
            "The shared INFERENCE_control_group value --crt-pool came from, recorded in the "
            "run's control-group provenance (see bin/control_group.py)"
        ),
    )
    parser.add_argument(
        "--control-group-provenance",
        default=None,
        choices=[None, "auto", "explicit", "per-method-override"],
        help="How --crt-pool was decided: from the declared MOI, from an explicit shared setting, or from INFERENCE_PERTURBO_CRT_POOL overriding it",
    )
    parser.add_argument(
        "--crt-test-control-elements",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Test the control elements too, so the control evaluation has p-values to score. "
            "Their own cells are in the pool they are tested against, which makes them conservative."
        ),
    )
    parser.add_argument("--device", default="gpu", help="JAX device for PerTurbo v2, e.g. gpu or cpu")
    parser.add_argument("--batch-size", type=int, default=0, help="SVI minibatch size")
    parser.add_argument("--num-steps-control", type=int, default=2500, help="Control-fit SVI steps")
    parser.add_argument(
        "--num-steps-betas", type=int, default=1000, help="Beta-fit SVI steps"
    )
    parser.add_argument(
        "--step-size",
        type=float,
        default=0.01,
        help=(
            "Adam learning rate for both SVI stages. 0.01 with 500 beta steps recovers simulated effects as "
            "well as 0.003 with 2,500 (PerTurbo's stage-two sweep, Sep 2026); 0.003 with 300 under-converges."
        ),
    )
    parser.add_argument("--max-chunk-size", type=int, default=50000, help="PerTurbo v2 max chunk cell count")
    parser.add_argument("--perturbation-chunk-size", type=int, default=0, help="PerTurbo v2 perturbation chunk size")
    parser.add_argument(
        "--crt-gene-chunk-size",
        type=int,
        default=0,
        help=(
            "Genes per conditional-randomization-test block; 0 lets the test size "
            "its own blocks. This is the test's memory knob and it is separate from "
            "--gene-chunk-size, which governs stage two. Our production Replogle "
            "runs use 500."
        ),
    )
    parser.add_argument(
        "--crt-max-gather-gib",
        type=float,
        default=0.0,
        help=(
            "Cap on one gather inside the test, in GiB; 0 leaves the default. "
            "Production Replogle runs use 8."
        ),
    )
    parser.add_argument(
        "--gene-chunk-size",
        type=int,
        default=0,
        help=(
            "Fit stage two in blocks of this many genes. 0 keeps the whole panel on "
            "device at once, which is what a large card can afford; a 40 GB card "
            "running a transcriptome-wide panel cannot, and asks for an allocation "
            "it will not get. Genes are independent given the baseline, so blocking "
            "changes the memory profile and not the result."
        ),
    )
    parser.add_argument("--size-factor-mode", default="observed", choices=["infer", "observed", "none"])
    parser.add_argument("--likelihood", default="negbin", choices=["nb", "negbin", "censored_nb", "lognormal_nb", "mixture_nb"])
    parser.add_argument("--prior", default="normal", choices=["normal", "cauchy"])
    parser.add_argument("--save-model-params", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--parallel-fits",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Fit the element and guide models concurrently (normally one GPU per fit).",
    )
    parser.add_argument("--element-gpu", default="0", help="CUDA device ID for the parallel element fit")
    parser.add_argument("--guide-gpu", default="1", help="CUDA device ID for the parallel guide fit")
    parser.add_argument(
        "--jax-cache-dir",
        default=".perturbo_jax_cache",
        help="Persistent JAX compilation cache root; separate element/guide subdirectories are used.",
    )
    return parser


def main() -> None:
    run_pipeline_adapter(build_parser().parse_args())


if __name__ == "__main__":
    main()
