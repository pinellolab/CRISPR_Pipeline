#!/usr/bin/env python
"""Build the SCEPTRE input for the non-targeting calibration pass.

PerTurbo tests every non-targeting pseudo-element (the ``non-targeting|N``
buckets) against every gene as part of its ordinary run, so its rows exist
already. SCEPTRE tests only the pairs it is handed and, under the non-targeting
contrast, refuses non-targeting targets outright. This step builds a second,
small SCEPTRE input in which those buckets are ordinary targets, so SCEPTRE
produces matched rows for the same hypotheses PerTurbo already has:

- the buckets are renamed ``ntc_calibration|N`` (guide ids are kept, so the
  rows join to PerTurbo's on ``(guide_id, gene_id)`` and on ``(N, gene_id)``);
- under ``complement`` nothing else changes;
- under ``nt_cells`` the NT-only cells are duplicated once more as the
  reference population, carrying one synthetic non-targeting gRNA, so each
  bucket is tested against a pool that contains its own cells, which is what
  PerTurbo's control-element test does under the control-anchored pool. Both
  tests are then mildly conservative by about n_k over N, and the summary
  records that ratio.

The pair sets: the main run's local pairs with each targeting element replaced
by a random bucket (``local``), every bucket crossed with a sample of all genes
(``global``), and every bucket crossed with the screen's direct-target genes
(``direct_target``), so the AUPRC check always has negatives on the positives'
genes. The main run's pairs table, families and outputs are not touched.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse

from analysis_output_formatting import make_h5mu_safe_dataframe
from inference_target_matching import normalize_target_identifier
from intended_target_key_utils import _control_mask, annotate_intended_target_groups

ASSIGNMENT_LAYER = "guide_assignment"
CALIBRATION_PREFIX = "ntc_calibration|"
CALIBRATION_TYPE = "ntc_calibration"
POOL_GUIDE_ID = "ntc_pool_copy"
POOL_CELL_SUFFIX = "__ntcpool"
PAIR_TYPE = "ntc_calibration"
PAIR_COLUMNS = [
    "guide_id", "gene_id", "intended_target_name", "intended_target_chr",
    "intended_target_start", "intended_target_end", "intended_target_key",
    "pair_type", "pair_set", "ntc_bucket",
]


# --- observations ------------------------------------------------------------


def _binary_assignment(guide) -> sparse.csr_matrix:
    matrix = guide.layers[ASSIGNMENT_LAYER] if ASSIGNMENT_LAYER in guide.layers else guide.X
    matrix = matrix if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    return (matrix > 0).astype(np.int8).tocsr()


def _bucket_of(name: str) -> int | None:
    text = str(name)
    if "|" not in text:
        return None
    try:
        return int(text.rsplit("|", 1)[1])
    except ValueError:
        return None


def guide_roles(guide_var: pd.DataFrame) -> pd.DataFrame:
    """One row per guide: its target, whether it is a control, and its NT bucket."""
    var = guide_var.copy()
    if "guide_id" not in var.columns:
        var["guide_id"] = var.index.astype(str)
    control = _control_mask(var).to_numpy(dtype=bool)
    names = var["intended_target_name"].astype(str)
    buckets = [(_bucket_of(name) if is_control else None) for name, is_control in zip(names, control)]
    return pd.DataFrame(
        {
            "guide_id": var["guide_id"].astype(str).to_numpy(),
            "intended_target_name": names.to_numpy(),
            "intended_target_key": var["intended_target_key"].astype(str).to_numpy(),
            "targeting": ~control,
            "role": np.where(control, "ntc", "targeting"),
            "ntc_bucket": pd.array(buckets, dtype="Int64"),
        }
    )


def cell_counts(guide, roles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Cells per guide and per element on the unrestricted data, with roles."""
    binary = _binary_assignment(guide)
    per_cell = np.asarray(binary.sum(axis=1)).ravel()
    control = (roles["role"] == "ntc").to_numpy()
    controls_per_cell = np.asarray(binary[:, control].sum(axis=1)).ravel()
    nt_carrying = controls_per_cell > 0
    nt_only = nt_carrying & (per_cell == controls_per_cell)

    n_cells = np.asarray(binary.sum(axis=0)).ravel()
    singleton_rows = sparse.diags((per_cell == 1).astype(np.int8)) @ binary
    n_singleton = np.asarray(singleton_rows.sum(axis=0)).ravel()
    per_guide = roles.copy()
    per_guide["n_cells"] = n_cells.astype(int)
    per_guide["n_cells_singleton"] = n_singleton.astype(int)

    keys = roles["intended_target_key"].to_numpy()
    codes, uniques = pd.factorize(keys, sort=False)
    indicator = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int8), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(uniques)),
    )
    per_element_cells = np.asarray(((binary @ indicator) > 0).sum(axis=0)).ravel()
    element_rows = (
        roles.groupby("intended_target_key", sort=False)
        .agg(element_name=("intended_target_name", "first"), role=("role", "first"),
             n_guides=("guide_id", "size"), ntc_bucket=("ntc_bucket", "first"))
        .reindex(uniques)
    )
    per_element = pd.DataFrame(
        {
            "intended_target_key": uniques,
            "element_name": element_rows["element_name"].to_numpy(),
            "role": element_rows["role"].to_numpy(),
            "n_guides": element_rows["n_guides"].to_numpy().astype(int),
            "ntc_bucket": element_rows["ntc_bucket"].to_numpy(),
            "n_cells": per_element_cells.astype(int),
        }
    )
    summary = {
        "n_cells": int(per_cell.size),
        "n_nt_carrying_cells": int(nt_carrying.sum()),
        "n_nt_only_cells": int(nt_only.sum()),
        "n_multi_guide_cells": int((per_cell > 1).sum()),
    }
    return per_guide, per_element, summary, nt_carrying, nt_only


# --- pair sets ----------------------------------------------------------------


def direct_target_pairs(roles: pd.DataFrame, gene_var: pd.DataFrame) -> pd.DataFrame:
    """(guide_id, gene_id) for every targeting guide and the gene it is aimed at,
    matched by gene id or by symbol: the positives of the AUPRC check."""
    ids = pd.Series(gene_var.index.astype(str), index=gene_var.index)
    keys = {"gene_id": normalize_target_identifier(ids)}
    for column in ("gene_name", "symbol", "gene_symbol"):
        if column in gene_var.columns:
            keys[column] = normalize_target_identifier(gene_var[column].astype("string"))
    lookup: dict[str, set[str]] = {}
    for column, values in keys.items():
        for gene_id, key in zip(ids, values):
            if pd.notna(key):
                lookup.setdefault(str(key), set()).add(str(gene_id))
    targeting = roles.loc[roles["targeting"], ["guide_id", "intended_target_name"]]
    targets = normalize_target_identifier(targeting["intended_target_name"].astype("string"))
    rows = []
    for guide_id, key in zip(targeting["guide_id"], targets):
        for gene_id in sorted(lookup.get(str(key), ())) if pd.notna(key) else ():
            rows.append((str(guide_id), gene_id))
    return pd.DataFrame(rows, columns=["guide_id", "gene_id"]).drop_duplicates().reset_index(drop=True)


def _direct_target_genes(roles: pd.DataFrame, gene_var: pd.DataFrame) -> list[str]:
    """Genes some targeting guide is aimed at, by id or by symbol."""
    return sorted(set(direct_target_pairs(roles, gene_var)["gene_id"]))


def build_element_pairs(
    local_pairs: pd.DataFrame,
    buckets: list[int],
    all_genes: list[str],
    direct_genes: list[str],
    *,
    max_local_pairs: int,
    max_global_pairs: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict]:
    """Element-level (bucket, gene) pairs with their set labels, and how they were sized."""
    sizing: dict = {}
    frames = []

    # local: the main run's targeting elements, each swapped for a random bucket.
    pairs = local_pairs.copy()
    pairs["intended_target_key"] = pairs["intended_target_key"].astype(str)
    pairs["gene_id"] = pairs["gene_id"].astype(str)
    is_control_pair = pairs["intended_target_key"].str.lower().str.contains("non-targeting", regex=False)
    pairs = pairs.loc[~is_control_pair, ["intended_target_key", "gene_id"]].drop_duplicates()
    gene_lists = pairs.groupby("intended_target_key", sort=False)["gene_id"].apply(list)
    sizing["local_elements_available"] = int(len(gene_lists))
    sizing["local_pairs_available"] = int(len(pairs))
    kept = list(gene_lists.index)
    if max_local_pairs > 0 and len(pairs) > max_local_pairs:
        order = rng.permutation(len(kept))
        total, selected = 0, []
        for position in order:
            element = kept[position]
            size = len(gene_lists[element])
            if total + size > max_local_pairs:
                continue
            selected.append(element)
            total += size
        kept = selected
    local_rows = []
    for element in kept:
        bucket = int(rng.choice(buckets))
        for gene in gene_lists[element]:
            local_rows.append((bucket, gene, element))
    local = pd.DataFrame(local_rows, columns=["ntc_bucket", "gene_id", "substituted_element"])
    local = local.drop_duplicates(subset=["ntc_bucket", "gene_id"])
    local["pair_set"] = "local"
    sizing["local_elements_used"] = int(len(kept))
    sizing["local_pairs"] = int(len(local))
    frames.append(local)

    # global: every bucket against a sample of all genes.
    if max_global_pairs > 0:
        n_genes = max(1, min(len(all_genes), max_global_pairs // max(1, len(buckets))))
    else:
        n_genes = len(all_genes)
    sampled = sorted(rng.choice(np.asarray(all_genes, dtype=object), size=n_genes, replace=False).tolist()) if n_genes < len(all_genes) else list(all_genes)
    global_pairs = pd.DataFrame(
        [(b, g, None) for b in buckets for g in sampled],
        columns=["ntc_bucket", "gene_id", "substituted_element"],
    )
    global_pairs["pair_set"] = "global"
    sizing["global_genes_sampled"] = int(len(sampled))
    sizing["global_pairs"] = int(len(global_pairs))
    frames.append(global_pairs)

    # direct_target: every bucket against the genes the targeting guides aim at.
    direct = pd.DataFrame(
        [(b, g, None) for b in buckets for g in direct_genes],
        columns=["ntc_bucket", "gene_id", "substituted_element"],
    )
    direct["pair_set"] = "direct_target"
    sizing["direct_target_genes"] = int(len(direct_genes))
    sizing["direct_target_pairs"] = int(len(direct))
    frames.append(direct)

    union = pd.concat(frames, ignore_index=True)
    union["ntc_bucket"] = union["ntc_bucket"].astype(int)
    labelled = (
        union.groupby(["ntc_bucket", "gene_id"], sort=False)
        .agg(
            pair_set=("pair_set", lambda s: "+".join(sorted(set(s)))),
            substituted_element=("substituted_element", lambda s: next((x for x in s if x is not None), None)),
        )
        .reset_index()
    )
    sizing["union_element_pairs"] = int(len(labelled))
    return labelled, sizing


# --- relabelling ---------------------------------------------------------------


def relabel_calibration_guides(guide_var: pd.DataFrame, roles: pd.DataFrame) -> pd.DataFrame:
    """Turn every NT bucket into an ordinary target named ntc_calibration|N."""
    var = guide_var.copy()
    if "guide_id" not in var.columns:
        var["guide_id"] = var.index.astype(str)
    is_ntc = (roles["role"] == "ntc").to_numpy()
    buckets = roles["ntc_bucket"].to_numpy()
    names = np.array([f"{CALIBRATION_PREFIX}{int(b)}" for b in buckets[is_ntc]], dtype=object)
    var = var.astype({c: object for c in ("intended_target_name", "intended_target_chr") if c in var.columns})
    var.loc[is_ntc, "intended_target_name"] = names
    # An empty chromosome, not a missing one: anndata cannot write a missing
    # value beside strings in an object column, and the key builder treats a
    # missing start or end as "no coordinates" regardless.
    var.loc[is_ntc, "intended_target_chr"] = ""
    var.loc[is_ntc, "intended_target_start"] = pd.NA
    var.loc[is_ntc, "intended_target_end"] = pd.NA
    var.loc[is_ntc, "targeting"] = True
    if "type" in var.columns:
        var["type"] = var["type"].astype(object)
        var.loc[is_ntc, "type"] = CALIBRATION_TYPE
    else:
        var["type"] = np.where(is_ntc, CALIBRATION_TYPE, "targeting")
    if "element_id" in var.columns:
        var["element_id"] = var["element_id"].astype(object)
        var.loc[is_ntc, "element_id"] = names
    return annotate_intended_target_groups(var)


def assert_relabel_invariants(var: pd.DataFrame, calibration_ids: set[str]) -> None:
    """What both methods' control detection must agree on, checked before writing."""
    from perturbo_v2_pipeline_adapter import _guide_control_mask

    ids = var["guide_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("guide ids are not unique after relabelling")
    is_cal = ids.isin(calibration_ids).to_numpy()
    pipeline_control = _control_mask(var).to_numpy(dtype=bool)
    adapter_control = _guide_control_mask(var).to_numpy(dtype=bool)
    if pipeline_control[is_cal].any() or adapter_control[is_cal].any():
        raise ValueError("a calibration pseudo-element still reads as a control")
    if not (pipeline_control[~is_cal].all() and adapter_control[~is_cal].all()):
        raise ValueError("the synthetic pool gRNA does not read as a control")
    for column in ("intended_target_name", "intended_target_key"):
        if var.loc[is_cal, column].astype(str).str.lower().str.contains("non-targeting", regex=False).any():
            raise ValueError(f"a calibration pseudo-element's {column} contains 'non-targeting'")
    if (var["intended_target_name"].astype(str) == "non-targeting").any():
        raise ValueError("a target is named exactly 'non-targeting'")


# --- assembly ------------------------------------------------------------------


def stored_obs_columns(path: str | Path) -> list[str]:
    """The top-level obs columns as written in the file, which is what the SCEPTRE
    driver reads as colData. ``mu.read_h5mu`` adds every modality column to the
    in-memory frame under a ``mod:`` prefix, and those must not be carried into
    an input whose every top-level column becomes a formula term."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        order = obs.attrs.get("column-order")
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in order] if order is not None else [
            k for k in obs.keys() if k != "_index"
        ]
    return [name for name in names if name != "_index"]


def _duplicate_rows(adata, rows: np.ndarray, suffix: str):
    copy = adata[rows].copy()
    copy.obs_names = [f"{name}{suffix}" for name in copy.obs_names]
    return copy


def assemble_mudata(full, roles: pd.DataFrame, relabelled_var: pd.DataFrame, pair_genes: list[str],
                    keep_cells: np.ndarray, pool_cells: np.ndarray, duplicate_pool: bool,
                    obs_columns: list[str]):
    """The SCEPTRE input: pair genes, calibration guides, kept cells, plus the pool copies."""
    gene = full["gene"]
    guide = full["guide"]
    gene_mask = gene.var_names.isin(pair_genes)
    is_ntc = (roles["role"] == "ntc").to_numpy()

    gene_kept = gene[keep_cells, gene_mask].copy()
    guide_kept = guide[keep_cells, is_ntc].copy()
    guide_kept.var = relabelled_var.loc[is_ntc].copy()
    guide_kept.var.index = guide_kept.var["guide_id"].astype(str).to_numpy()
    obs_kept = full.obs.loc[gene_kept.obs_names, obs_columns].copy()

    if duplicate_pool:
        copy_rows = keep_cells & pool_cells
        gene_copy = _duplicate_rows(gene, copy_rows, POOL_CELL_SUFFIX)
        gene_copy = gene_copy[:, gene_mask].copy()
        n_copy = gene_copy.n_obs
        pool_var = pd.DataFrame(
            {
                "guide_id": [POOL_GUIDE_ID],
                "targeting": [False],
                "type": ["non-targeting"],
                "intended_target_name": ["non-targeting|pool"],
                "intended_target_chr": [""],
                "intended_target_start": [pd.NA],
                "intended_target_end": [pd.NA],
            },
            index=[POOL_GUIDE_ID],
        )
        # Calibration guides: zero on the copies. Pool gRNA: one on every copy.
        zeros = sparse.csr_matrix((n_copy, guide_kept.n_vars), dtype=np.float32)
        ones = sparse.csr_matrix(np.ones((n_copy, 1), dtype=np.float32))
        guide_copy = ad.AnnData(
            X=sparse.hstack([zeros, ones]).tocsr(),
            obs=gene_copy.obs[[]].copy(),
            var=pd.concat([guide_kept.var[["guide_id", "targeting", "type", "intended_target_name",
                                           "intended_target_chr", "intended_target_start", "intended_target_end"]],
                           pool_var]),
        )
        guide_copy.layers[ASSIGNMENT_LAYER] = guide_copy.X.copy()
        guide_var = pd.concat([guide_kept.var, pool_var])
        guide_var = annotate_intended_target_groups(guide_var)
        guide_var.index = guide_var["guide_id"].astype(str).to_numpy()
        # Originals: zero on the pool gRNA.
        orig_assignment = _binary_assignment(guide_kept).astype(np.float32)
        orig_counts = guide_kept.X if sparse.issparse(guide_kept.X) else sparse.csr_matrix(guide_kept.X)
        pool_zero = sparse.csr_matrix((guide_kept.n_obs, 1), dtype=np.float32)
        guide_all = ad.AnnData(
            X=sparse.vstack([sparse.hstack([orig_counts.astype(np.float32), pool_zero]), guide_copy.X]).tocsr(),
            obs=pd.concat([guide_kept.obs[[]], guide_copy.obs[[]]]),
            var=guide_var,
        )
        guide_all.layers[ASSIGNMENT_LAYER] = sparse.vstack(
            [sparse.hstack([orig_assignment, pool_zero]), guide_copy.layers[ASSIGNMENT_LAYER]]
        ).tocsr()
        for key, value in guide.uns.items():
            guide_all.uns[key] = value
        gene_all = ad.concat([gene_kept, gene_copy], join="outer", merge="same")
        gene_all.uns = dict(gene_kept.uns)
        obs_copy = obs_kept.loc[gene_kept.obs_names[copy_rows[keep_cells]]].copy()
        obs_copy.index = [f"{name}{POOL_CELL_SUFFIX}" for name in obs_copy.index]
        obs_all = pd.concat([obs_kept, obs_copy])
    else:
        guide_var = annotate_intended_target_groups(guide_kept.var)
        guide_var.index = guide_var["guide_id"].astype(str).to_numpy()
        guide_all = guide_kept
        guide_all.var = guide_var
        if ASSIGNMENT_LAYER not in guide_all.layers:
            guide_all.layers[ASSIGNMENT_LAYER] = _binary_assignment(guide_all).astype(np.float32)
        gene_all = gene_kept
        obs_all = obs_kept

    out = mu.MuData({"gene": gene_all, "guide": guide_all})
    out.uns = dict(full.uns)
    out.obs = obs_all.loc[out.obs_names].copy()
    return out


# --- main ------------------------------------------------------------------------


def _write_skip(args, reason: str, summary: dict, per_guide: pd.DataFrame | None, per_element: pd.DataFrame | None) -> int:
    summary = {**summary, "status": "skipped", "reason": reason}
    Path(args.summary).write_text(json.dumps(summary, indent=2, default=str) + "\n")
    if per_guide is not None:
        per_guide.to_csv(args.cell_counts_per_guide, sep="\t", index=False)
        per_element.to_csv(args.cell_counts_per_element, sep="\t", index=False)
        per_guide[["guide_id", "intended_target_name", "intended_target_key", "targeting", "role", "ntc_bucket"]].to_csv(args.guide_targets, sep="\t", index=False)
    pd.DataFrame(columns=["guide_id", "gene_id"]).to_csv(args.direct_targets, sep="\t", index=False)
    print(f"NTC calibration skipped: {reason}", flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rng = np.random.default_rng(int(args.seed))
    regime = args.sceptre_control_group
    if regime not in ("nt_cells", "complement"):
        raise SystemExit(f"--sceptre-control-group must be nt_cells or complement, got {regime!r}")
    duplicate_pool = regime == "nt_cells"

    print(f"Reading {args.full_mudata}...", flush=True)
    full = mu.read_h5mu(args.full_mudata)
    full["guide"].var = annotate_intended_target_groups(full["guide"].var)
    roles = guide_roles(full["guide"].var)
    base_summary = {
        "regime": regime,
        "analysis_moi": args.analysis_moi,
        "seed": int(args.seed),
        "max_local_pairs": int(args.max_local_pairs),
        "max_global_pairs": int(args.max_global_pairs),
        "n_guides": int(len(roles)),
        "n_ntc_guides": int((roles["role"] == "ntc").sum()),
    }
    per_guide, per_element, counts, nt_carrying, nt_only = cell_counts(full["guide"], roles)
    base_summary.update(counts)
    buckets = sorted(int(b) for b in roles.loc[roles["role"] == "ntc", "ntc_bucket"].dropna().unique())
    base_summary["n_ntc_buckets"] = len(buckets)
    if not buckets:
        return _write_skip(args, "the screen has no non-targeting guides", base_summary, per_guide, per_element)
    if duplicate_pool and counts["n_nt_only_cells"] == 0:
        return _write_skip(args, "no cell carries only non-targeting guides, so there is no pool to duplicate", base_summary, per_guide, per_element)

    print(f"Reading the local pairs from {args.pairs_mudata}...", flush=True)
    pairs_source = mu.read_h5mu(args.pairs_mudata, backed="r")
    local_pairs = pairs_source.uns.get("pairs_to_test")
    if local_pairs is None:
        raise SystemExit("the pairs MuData carries no uns['pairs_to_test']")
    local_pairs = local_pairs.copy() if isinstance(local_pairs, pd.DataFrame) else pd.DataFrame(local_pairs)
    if "intended_target_key" not in local_pairs.columns:
        key_by_guide = dict(zip(roles["guide_id"], roles["intended_target_key"]))
        local_pairs["intended_target_key"] = local_pairs["guide_id"].astype(str).map(key_by_guide)
    file_manager = getattr(pairs_source, "file", None)
    if file_manager is not None:
        file_manager.close()

    all_genes = [str(g) for g in full["gene"].var_names]
    direct = direct_target_pairs(roles, full["gene"].var)
    direct_genes = sorted(set(direct["gene_id"]))
    element_pairs, sizing = build_element_pairs(
        local_pairs, buckets, all_genes, direct_genes,
        max_local_pairs=int(args.max_local_pairs), max_global_pairs=int(args.max_global_pairs), rng=rng,
    )
    base_summary.update(sizing)

    relabelled = relabel_calibration_guides(full["guide"].var, roles)
    calibration_ids = set(roles.loc[roles["role"] == "ntc", "guide_id"])
    # Guide-level expansion over each bucket's guides, with the relabelled metadata.
    meta = relabelled.loc[relabelled["guide_id"].astype(str).isin(calibration_ids)].copy()
    meta["ntc_bucket"] = meta["intended_target_name"].astype(str).map(_bucket_of).astype(int)
    guide_pairs = element_pairs.merge(
        meta[["guide_id", "ntc_bucket", "intended_target_name", "intended_target_chr",
              "intended_target_start", "intended_target_end", "intended_target_key"]],
        on="ntc_bucket", how="inner",
    )
    guide_pairs["pair_type"] = PAIR_TYPE
    guide_pairs = guide_pairs[PAIR_COLUMNS + ["substituted_element"]].sort_values(["ntc_bucket", "gene_id", "guide_id"]).reset_index(drop=True)
    base_summary["union_guide_pairs"] = int(len(guide_pairs))

    keep_cells = nt_carrying if duplicate_pool else np.ones(full.n_obs, dtype=bool)
    pair_genes = sorted(set(guide_pairs["gene_id"]))
    obs_columns = stored_obs_columns(args.full_mudata)
    out = assemble_mudata(full, roles, relabelled, pair_genes, keep_cells, nt_only, duplicate_pool, obs_columns)
    assert_relabel_invariants(out["guide"].var, calibration_ids)
    if list(out.obs.columns) != obs_columns:
        raise ValueError(f"top-level obs columns changed: {list(out.obs.columns)} vs {obs_columns}")
    base_summary["sceptre_covariate_columns"] = obs_columns
    out.uns["pairs_to_test"] = make_h5mu_safe_dataframe(guide_pairs[PAIR_COLUMNS])

    # Leniency: how much of its own reference each bucket is.
    pool_size = int(counts["n_nt_only_cells"]) if duplicate_pool else int(full.n_obs)
    bucket_cells = per_element.loc[per_element["role"] == "ntc", ["ntc_bucket", "n_cells"]].dropna()
    base_summary["self_inclusion_ratio_by_bucket"] = {
        str(int(b)): (float(n) / pool_size if pool_size else None) for b, n in zip(bucket_cells["ntc_bucket"], bucket_cells["n_cells"])
    }
    base_summary.update({
        "status": "ok",
        "n_cells_kept": int(out["gene"].n_obs if not duplicate_pool else keep_cells.sum()),
        "n_pool_copies": int((keep_cells & nt_only).sum()) if duplicate_pool else 0,
        "n_cells_in_input": int(out.n_obs),
        "n_genes_in_input": int(out["gene"].n_vars),
        "n_guides_in_input": int(out["guide"].n_vars),
        "reference_pool_size": pool_size,
    })

    print(f"Writing {args.output_mudata} ({out.n_obs:,} cells, {out['gene'].n_vars:,} genes, {out['guide'].n_vars:,} guides)...", flush=True)
    out.write(args.output_mudata)
    guide_pairs.to_csv(args.pairs, sep="\t", index=False, compression="gzip")
    per_guide.to_csv(args.cell_counts_per_guide, sep="\t", index=False)
    per_element.to_csv(args.cell_counts_per_element, sep="\t", index=False)
    per_guide[["guide_id", "intended_target_name", "intended_target_key", "targeting", "role", "ntc_bucket"]].to_csv(args.guide_targets, sep="\t", index=False)
    direct.to_csv(args.direct_targets, sep="\t", index=False)
    base_summary["direct_target_guide_pairs"] = int(len(direct))
    Path(args.summary).write_text(json.dumps(base_summary, indent=2, default=str) + "\n")
    print(
        f"NTC calibration input ready: {len(buckets)} bucket(s), {sizing['union_element_pairs']:,} element pairs "
        f"({sizing['local_pairs']:,} local, {sizing['global_pairs']:,} global, {sizing['direct_target_pairs']:,} direct-target), "
        f"{len(guide_pairs):,} guide pairs, regime {regime}.",
        flush=True,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--full-mudata", required=True, help="prepared MuData with every guide (the one PerTurbo fits)")
    parser.add_argument("--pairs-mudata", required=True, help="prepared MuData carrying uns['pairs_to_test'] (the local pairs)")
    parser.add_argument("--sceptre-control-group", required=True, choices=["nt_cells", "complement"])
    parser.add_argument("--analysis-moi", default="", help="recorded only; the SCEPTRE process receives it separately")
    parser.add_argument("--max-local-pairs", type=int, default=10000)
    parser.add_argument("--max-global-pairs", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-mudata", default="ntc_calibration_input.h5mu")
    parser.add_argument("--pairs", default="ntc_calibration_pairs.tsv.gz")
    parser.add_argument("--guide-targets", default="ntc_guide_targets.tsv")
    parser.add_argument("--direct-targets", default="ntc_direct_targets.tsv", help="(guide_id, gene_id) positives for the AUPRC check")
    parser.add_argument("--cell-counts-per-guide", default="ntc_cell_counts_per_guide.tsv")
    parser.add_argument("--cell-counts-per-element", default="ntc_cell_counts_per_element.tsv")
    parser.add_argument("--summary", default="ntc_calibration_design.json")
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
