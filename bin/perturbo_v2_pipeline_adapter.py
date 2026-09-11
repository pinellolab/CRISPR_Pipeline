#!/usr/bin/env python
"""Run PerTurbo v2 while preserving CRISPR_Pipeline's inference outputs."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control

from analysis_output_formatting import make_h5mu_safe_dataframe
from intended_target_key_utils import (
    annotate_intended_target_groups,
    enrich_pairs_with_target_metadata,
    get_target_lookup,
)
from mudata_uns_io import write_uns_patch
from result_table_io import write_result_table


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
    guide_obs = mdata[GUIDE_MODALITY].obs
    gene_obs = mdata[GENE_MODALITY].obs
    if "total_guide_umis" in guide_obs.columns:
        values = np.asarray(guide_obs["total_guide_umis"], dtype=float)
    else:
        values = np.asarray(_get_assignment_matrix(mdata).sum(axis=1)).ravel()
        guide_obs["total_guide_umis"] = values
    centered = np.log1p(values)
    centered = centered - float(np.nanmean(centered))
    gene_obs["log1p_total_guide_umis_centered"] = centered


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


def _covariate_has_control_variance(input_path: Path, map_key: str, names_key: str) -> bool:
    with _open_mudata(input_path, backed="r") as mdata:
        values = pd.to_numeric(
            mdata[GENE_MODALITY].obs["log1p_total_guide_umis_centered"],
            errors="coerce",
        ).to_numpy(dtype=float)
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

    # PerTurbo uses fitted control-element z-values to calibrate empirical
    # p-values. Include every non-targeting element for each gene in the
    # requested hypothesis set, without expanding back to all RNA genes.
    control_elements = [name for name in element_names if CONTROL_SUBSTRING in name.lower()]
    tested_genes = requested["gene"].drop_duplicates().tolist()
    if control_elements:
        controls = pd.MultiIndex.from_product(
            [control_elements, tested_genes], names=["element", "gene"]
        ).to_frame(index=False)
        native_pairs = pd.concat([requested, controls], ignore_index=True).drop_duplicates(ignore_index=True)
    else:
        native_pairs = requested
    print(
        f"Prepared native PerTurbo pairs-to-test table: {len(requested):,} requested pairs plus "
        f"{len(control_elements) * len(tested_genes):,} control-null pairs "
        f"({len(native_pairs):,} unique total)."
    )
    return native_pairs


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
) -> pd.DataFrame:
    with _open_mudata(prepared_mudata_path, backed="r") as mdata:
        target_lookup = get_target_lookup(pd.DataFrame(mdata[GUIDE_MODALITY].var).copy())
    out = _convert_common_effect_columns(effects).rename(columns={"element": "intended_target_key"})
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
    ]
    out["perturbo_q_value"] = _bh_adjust(out["p_value"])
    return out


def convert_guide_effects(
    effects: pd.DataFrame,
    guide_name_map: dict[str, str],
    prepared_mudata_path: str | Path,
    *,
    test_all_pairs: bool,
) -> pd.DataFrame:
    out = _convert_common_effect_columns(effects)
    out["guide_id"] = out["element"].map(guide_name_map).fillna(out["element"]).astype(str)
    out = _maybe_filter_pairs(out, prepared_mudata_path, inference_type="guide", test_all_pairs=test_all_pairs)
    out = out[["gene_id", "guide_id", "log2_fc", "perturbo_fc_se", "p_value"]]
    out["perturbo_q_value"] = _bh_adjust(out["p_value"])
    return out


def _convert_common_effect_columns(effects: pd.DataFrame) -> pd.DataFrame:
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
    out["p_value"] = p_value.where(p_value.notna(), posterior)
    out["perturbo_posterior_prob"] = posterior
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
    if _covariate_has_control_variance(input_path, map_key, names_key):
        cmd.extend(["--continuous-covariates", "log1p_total_guide_umis_centered"])
    else:
        print(
            "Skipping log1p_total_guide_umis_centered covariate because it has "
            "zero variance in PerTurbo control cells for this run."
        )
    with _open_mudata(input_path, backed="r") as mdata:
        has_batch = "batch" in mdata[GENE_MODALITY].obs.columns
    if has_batch:
        cmd.extend(["--batch-covariate", "batch"])
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


def run_pipeline_adapter(args: argparse.Namespace) -> None:
    with tempfile.TemporaryDirectory(prefix="perturbo_v2_pipeline_") as tmp:
        tmp_dir = Path(tmp)
        prepared = tmp_dir / "prepared_mudata.h5mu"
        guide_name_map = prepare_mudata_for_perturbo_v2(
            args.input,
            prepared,
            test_all_pairs=args.test_all_pairs,
        )

        args.resolved_crt_pool = _resolve_crt_pool(args.crt_pool, _read_moi(args.input, args.moi))
        if args.crt:
            print(f"PerTurbo CRT pool: {args.resolved_crt_pool} (requested {args.crt_pool}).")
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

        element_dir = tmp_dir / "element_fit"
        guide_dir = tmp_dir / "guide_fit"
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
        if not args.parallel_fits:
            _wait_for_fits({"element": element_process})
        guide_process = _run_perturbo(
            prepared,
            guide_dir,
            map_key=GUIDE_MAP_KEY,
            names_key=GUIDE_NAMES_KEY,
            pairs_to_test_path=guide_pairs_path,
            gpu_id=args.guide_gpu if args.parallel_fits else None,
            phase="guide",
            args=args,
        )
        _wait_for_fits({"element": element_process, "guide": guide_process} if args.parallel_fits else {"guide": guide_process})

        element_effects = pd.read_parquet(element_dir / "element_effects.parquet")
        guide_effects = pd.read_parquet(guide_dir / "element_effects.parquet")
        # Every pair, one Benjamini-Hochberg family: the transcriptome-wide tables.
        global_element_df = convert_element_effects(element_effects, prepared, test_all_pairs=True)
        global_guide_df = convert_guide_effects(guide_effects, guide_name_map, prepared, test_all_pairs=True)
        # The requested pairs, corrected within that set alone: the local tables.
        local_element_df = local_guide_df = None
        if wants_local:
            local_element_df = convert_element_effects(
                _requested_pairs_table(element_dir, element_effects, element_pairs_path), prepared, test_all_pairs=True
            )
            local_guide_df = convert_guide_effects(
                _requested_pairs_table(guide_dir, guide_effects, guide_pairs_path), guide_name_map, prepared, test_all_pairs=True
            )
        # --per-*-output keep their historical meaning: global with --test-all-pairs,
        # local without. The --local-per-*-output files add the local tables beside
        # the global ones so a single run serves both.
        primary_element_df = global_element_df if args.test_all_pairs else local_element_df
        primary_guide_df = global_guide_df if args.test_all_pairs else local_guide_df
        write_result_table(primary_element_df, args.per_element_output)
        write_result_table(primary_guide_df, args.per_guide_output)
        if args.local_per_element_output:
            write_result_table(local_element_df, args.local_per_element_output)
        if args.local_per_guide_output:
            write_result_table(local_guide_df, args.local_per_guide_output)

        if args.output_mudata:
            updates = {
                # Keep the generic keys used by standalone/single-method consumers,
                # and expose the analysis-qualified keys consumed by the pipeline's
                # evaluation and dashboard steps.
                "per_element_results": make_h5mu_safe_dataframe(primary_element_df),
                "per_guide_results": make_h5mu_safe_dataframe(primary_guide_df),
                "global_analysis_per_element_results": make_h5mu_safe_dataframe(global_element_df),
                "global_analysis_per_guide_results": make_h5mu_safe_dataframe(global_guide_df),
            }
            if local_element_df is not None:
                updates["local_analysis_per_element_results"] = make_h5mu_safe_dataframe(local_element_df)
                updates["local_analysis_per_guide_results"] = make_h5mu_safe_dataframe(local_guide_df)
            write_uns_patch(args.input, args.output_mudata, updates=updates)

        if args.v2_artifact_dir:
            artifact_dir = Path(args.v2_artifact_dir)
            artifact_dir.mkdir(parents=True, exist_ok=True)
            for name, source in {
                "element_pairs_to_test.parquet": element_pairs_path,
                "guide_pairs_to_test.parquet": guide_pairs_path,
            }.items():
                if source is not None:
                    shutil.copy2(source, artifact_dir / name)
            for name, source in {"element": element_dir, "guide": guide_dir}.items():
                target = artifact_dir / name
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run PerTurbo v2 and emit CRISPR_Pipeline-compatible result files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input CRISPR_Pipeline MuData file")
    parser.add_argument("--per-element-output", required=True, help="Output per-element table (.tsv.gz or .parquet, by extension)")
    parser.add_argument("--per-guide-output", required=True, help="Output per-guide table (.tsv.gz or .parquet, by extension)")
    parser.add_argument("--output-mudata", help="Optional MuData with result tables in .uns")
    parser.add_argument("--v2-artifact-dir", default=None, help="Optional directory for raw PerTurbo v2 artifacts")
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
    parser.add_argument("--crt", action=argparse.BooleanOptionalAction, default=True, help="Run the conditional randomization test (exact saddlepoint, no resampling)")
    parser.add_argument(
        "--crt-pool",
        default="from-moi",
        choices=["from-moi", "auto", "all-cells", "control-anchored"],
        help="Cells a perturbation is tested against; from-moi maps the pipeline's MOI setting (high -> all-cells, low -> control-anchored)",
    )
    parser.add_argument("--moi", default=None, help="Override the MOI setting read from guide.uns['moi']")
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
