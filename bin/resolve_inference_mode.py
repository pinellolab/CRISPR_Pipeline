#!/usr/bin/env python
"""Settle the inference mode once, from the data, for SCEPTRE and PerTurbo alike.

Reads the MuData PerTurbo fits (every guide, the analysed cells), measures the
guide assignments, and resolves the shared control-group setting into explicit
values for both methods plus the MOI mode the SCEPTRE driver must use. The
result is one JSON file the inference subworkflow hands to both processes, so
the two methods can no longer settle the design differently on their own: the
PerTurbo adapter used to override a declared MOI from the assignments it saw,
and the SCEPTRE driver applied the same rule to the cis subset without ever
changing its control group.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mudata as md
import numpy as np
from scipy import sparse

from control_group import ControlGroupError, normalize_moi, resolve_from_observations
from intended_target_key_utils import _control_mask

ASSIGNMENT_LAYER = "guide_assignment"


def _declared_moi_from_uns(guide) -> str | None:
    raw = guide.uns.get("moi") if hasattr(guide, "uns") else None
    if raw is None:
        return None
    values = np.asarray(raw).ravel()
    if values.size == 0:
        return None
    first = values[0]
    return first.decode() if isinstance(first, (bytes, np.bytes_)) else str(first)


def observe_assignments(mdata) -> dict:
    """Per-cell guide counts and the control-only population, on the full matrix."""
    guide = mdata["guide"]
    matrix = guide.layers[ASSIGNMENT_LAYER] if ASSIGNMENT_LAYER in guide.layers else guide.X
    matrix = matrix if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    binary = (matrix > 0).tocsr()
    per_cell = np.asarray(binary.sum(axis=1)).ravel()
    control = _control_mask(guide.var).to_numpy(dtype=bool)
    controls_per_cell = np.asarray(binary[:, control].sum(axis=1)).ravel()
    targeting_per_cell = per_cell - controls_per_cell
    control_only = (controls_per_cell > 0) & (targeting_per_cell == 0)
    return {
        "n_cells": int(per_cell.size),
        "n_guides": int(binary.shape[1]),
        "n_control_guides": int(control.sum()),
        "max_guides_per_cell": int(per_cell.max()) if per_cell.size else 0,
        "n_multi_guide_cells": int((per_cell > 1).sum()),
        "n_unassigned_cells": int((per_cell == 0).sum()),
        "n_control_only_cells": int(control_only.sum()),
        "observed_singleton": bool(per_cell.size and per_cell.max() <= 1),
        "declared_moi_in_mudata": _declared_moi_from_uns(guide),
    }


def resolve(mudata_path: str | Path, *, setting, declared_moi, perturbo_pool, sceptre_group) -> dict:
    mdata = md.read_h5mu(mudata_path, backed="r")
    try:
        observed = observe_assignments(mdata)
    finally:
        file_manager = getattr(mdata, "file", None)
        if file_manager is not None:
            file_manager.close()
    declared = declared_moi if declared_moi not in (None, "") else observed["declared_moi_in_mudata"]
    resolution = resolve_from_observations(
        setting,
        declared,
        observed_singleton=observed["observed_singleton"],
        n_control_only_cells=observed["n_control_only_cells"],
        n_multi_guide_cells=observed["n_multi_guide_cells"],
        perturbo_pool=perturbo_pool or None,
        sceptre_group=sceptre_group or None,
    )
    resolution["declared_moi_parameter"] = normalize_moi(declared_moi)
    resolution["observations"] = observed
    resolution["mudata"] = str(mudata_path)
    return resolution


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mudata", required=True, help="MuData with every guide and the analysed cells")
    parser.add_argument("--setting", default="auto", help="INFERENCE_control_group")
    parser.add_argument("--declared-moi", default="", help="Multiplicity_of_infection; blank reads guide.uns['moi']")
    parser.add_argument("--perturbo-pool", default="", help="INFERENCE_PERTURBO_CRT_POOL (per-method override)")
    parser.add_argument("--sceptre-group", default="", help="INFERENCE_SCEPTRE_control_group (per-method override)")
    parser.add_argument("--output", default="inference_mode.json")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        resolution = resolve(
            args.mudata,
            setting=args.setting,
            declared_moi=args.declared_moi,
            perturbo_pool=args.perturbo_pool,
            sceptre_group=args.sceptre_group,
        )
    except ControlGroupError as error:
        raise SystemExit(f"ERROR: {error}") from error
    Path(args.output).write_text(json.dumps(resolution, indent=2, default=str) + "\n")
    print(resolution["log_line"], flush=True)
    for note in resolution["notes"]:
        print(f"Note: {note}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
