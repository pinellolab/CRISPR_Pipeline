#!/usr/bin/env python3
"""Validate a pre-inference handoff and apply external cell/covariate edits."""

import argparse
import hashlib
import json
import re
from pathlib import Path

import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype={"cell_id": str})


def matrix_equal(left, right):
    if left.shape != right.shape:
        return False
    if sparse.issparse(left) or sparse.issparse(right):
        delta = sparse.csr_matrix(left) - sparse.csr_matrix(right)
        return delta.nnz == 0 or np.all(delta.data == 0)
    return np.array_equal(np.asarray(left), np.asarray(right))


def validate_edited(parent, edited):
    if not {"gene", "guide"}.issubset(edited.mod):
        raise ValueError("edited MuData requires gene and guide modalities")
    parent_cells = parent["gene"].obs_names
    cells = edited["gene"].obs_names
    if not cells.is_unique or not edited["guide"].obs_names.equals(cells):
        raise ValueError("edited gene/guide cell identifiers are not identical and unique")
    if not cells.isin(parent_cells).all():
        raise ValueError("edited MuData adds cells absent from the parent checkpoint")
    for modality in ("gene", "guide"):
        if not edited[modality].var_names.equals(parent[modality].var_names):
            raise ValueError(f"edited {modality} feature IDs differ from the parent")
        original = parent[modality][cells, :]
        if not matrix_equal(original.X, edited[modality].X):
            raise ValueError(f"edited {modality}.X differs from parent raw counts")
    parent_layers = parent["guide"].layers
    edited_layers = edited["guide"].layers
    if "guide_assignment" not in parent_layers or "guide_assignment" not in edited_layers:
        raise ValueError("guide_assignment layer is required in both MuData files")
    if not matrix_equal(parent["guide"][cells, :].layers["guide_assignment"], edited_layers["guide_assignment"]):
        raise ValueError("edited guide_assignment differs from the parent")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-manifest", required=True, type=Path)
    parser.add_argument("--parent-mudata", required=True, type=Path)
    parser.add_argument("--input-mudata", type=Path)
    parser.add_argument("--keep-cells", type=Path)
    parser.add_argument("--covariates", type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", default="analysis_input.h5mu", type=Path)
    parser.add_argument("--report", default="external_mudata_report.json", type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.parent_manifest.read_text())
    if manifest.get("status") != "READY_FOR_INFERENCE":
        raise ValueError("parent manifest is not a ready pre-inference checkpoint")
    if sha256(args.parent_mudata) != manifest.get("mudata_sha256"):
        raise ValueError("parent checkpoint SHA-256 does not match the manifest")
    if sha256(args.gtf) != manifest.get("gtf_sha256"):
        raise ValueError("reference GTF differs from the checkpoint run")

    parent = mu.read_h5mu(args.parent_mudata)
    source = args.input_mudata or args.parent_mudata
    edited = mu.read_h5mu(source)
    validate_edited(parent, edited)
    starting_cells = edited["gene"].obs_names.copy()
    external_covariates = []

    if args.keep_cells:
        table = read_table(args.keep_cells)
        if list(table.columns) != ["cell_id"]:
            raise ValueError("cell allowlist must have exactly one column named cell_id")
        requested = pd.Index(table["cell_id"].astype(str))
        if not requested.is_unique or requested.isna().any():
            raise ValueError("cell allowlist contains duplicate or missing IDs")
        unknown = requested.difference(starting_cells)
        if len(unknown):
            raise ValueError(f"cell allowlist contains {len(unknown)} IDs absent from the input MuData")
        if len(requested) == 0:
            raise ValueError("cell allowlist is empty")
        retained = starting_cells[starting_cells.isin(requested)]
        edited = edited[retained, :].copy()

    if args.covariates:
        table = read_table(args.covariates)
        if "cell_id" not in table or len(table.columns) < 2:
            raise ValueError("covariate table needs cell_id and at least one covariate")
        ids = pd.Index(table["cell_id"].astype(str))
        if not ids.is_unique or ids.isna().any():
            raise ValueError("covariate table contains duplicate or missing cell IDs")
        cells = edited["gene"].obs_names
        if not cells.isin(ids).all():
            raise ValueError("covariate table is missing retained cells")
        table.index = ids
        for name in table.columns.drop("cell_id"):
            if not re.fullmatch(r"external_[A-Za-z0-9_]+", name):
                raise ValueError(f"custom covariate {name!r} must begin with external_")
            if any(name in mod.obs for mod in edited.mod.values()) or name in edited.obs:
                raise ValueError(f"custom covariate {name!r} would overwrite an existing column")
            values = pd.to_numeric(table.loc[cells, name], errors="raise").to_numpy(dtype=float)
            if not np.isfinite(values).all() or np.unique(values).size < 2:
                raise ValueError(f"custom covariate {name!r} must be finite and nonconstant")
            edited["gene"].obs[name] = values
            edited.obs[name] = values
            external_covariates.append(name)

    if edited["gene"].n_obs == 0:
        raise ValueError("external selection removed every cell")
    edited.uns["external_covariates"] = external_covariates
    edited.uns["preinference_parent_sha256"] = manifest["mudata_sha256"]
    edited.write_h5mu(args.output)
    report = {
        "parent_run_name": manifest["run_name"],
        "parent_sha256": manifest["mudata_sha256"],
        "input_sha256": sha256(source),
        "output_sha256": sha256(args.output),
        "parent_cells": int(parent["gene"].n_obs),
        "submitted_cells": int(len(starting_cells)),
        "accepted_cells": int(edited["gene"].n_obs),
        "external_covariates": external_covariates,
        "filter_file": args.keep_cells.name if args.keep_cells else None,
        "covariate_file": args.covariates.name if args.covariates else None,
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
