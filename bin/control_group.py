#!/usr/bin/env python
"""One setting for the cells a perturbation is compared against.

SCEPTRE and PerTurbo make the same statistical choice under different names.
``INFERENCE_control_group`` states it once, in SCEPTRE's vocabulary, and this
module is the single place the mapping lives:

    setting        SCEPTRE control_group   PerTurbo --crt-pool
    -----------    ---------------------   -------------------
    nt_cells       nt_cells                control-anchored
    complement     complement              all-cells
    auto           from the declared MOI   from the declared MOI

``auto`` reads ``Multiplicity_of_infection``: ``low`` means the non-targeting
cells are the reference population, ``high`` means the complement of the
perturbation is (the only contrast SCEPTRE offers at high MOI). When the MOI is
neither, each method measures the design from the data as far as it can --
PerTurbo does that natively (``--crt-pool auto``), SCEPTRE cannot, so it falls
back to ``complement``, which is valid for every design.

The Groovy mirror of this table lives in ``modules/local/control_group/main.nf``
and is what the pipeline actually passes to the two processes;
``tests/test_control_group_resolution.py`` asserts the two agree.
"""

from __future__ import annotations

import argparse
import json

SETTING_PARAM = "INFERENCE_control_group"
PERTURBO_PARAM = "INFERENCE_PERTURBO_CRT_POOL"
SCEPTRE_PARAM = "INFERENCE_SCEPTRE_control_group"

#: The values ``INFERENCE_control_group`` accepts.
SETTINGS = ("auto", "nt_cells", "complement")

#: The historical per-method defaults. A per-method parameter still sitting on
#: its historical value is "not set" for precedence purposes: every config
#: written before the shared setting existed carries these.
PERTURBO_HISTORICAL_DEFAULT = "from-moi"
SCEPTRE_HISTORICAL_DEFAULT = "complement"

#: What each setting means to each method. ``from-moi`` is PerTurbo's own name
#: for "map the declared MOI", and keeping it means the adapter's observed-design
#: override still applies, exactly as it did before this setting existed.
SCEPTRE_BY_SETTING = {"nt_cells": "nt_cells", "complement": "complement"}
PERTURBO_BY_SETTING = {"nt_cells": "control-anchored", "complement": "all-cells"}

#: Where ``from-moi`` lands, for the log line and the recorded provenance. The
#: adapter computes this itself from the same table; this is only a report.
PERTURBO_BY_MOI = {"low": "control-anchored", "high": "all-cells"}
SCEPTRE_BY_MOI = {"low": "nt_cells", "high": "complement"}


class ControlGroupError(ValueError):
    """An unusable control-group configuration."""


def normalize_moi(moi: str | None) -> str:
    """The declared MOI, or ``'unknown'`` when it is absent or something else."""
    text = "" if moi is None else str(moi).strip().lower()
    return text if text in ("low", "high") else "unknown"


def normalize_setting(setting: str | None) -> str:
    text = "auto" if setting is None else str(setting).strip().lower()
    if text not in SETTINGS:
        raise ControlGroupError(
            f"{SETTING_PARAM}={setting!r} is not a control group. "
            f"Accepted values: {', '.join(SETTINGS)}. "
            "'nt_cells' compares each perturbation with the non-targeting cells, "
            "'complement' with every other cell, and 'auto' takes whichever the "
            "declared Multiplicity_of_infection implies."
        )
    return text


def _is_override(value: str | None, historical_default: str) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return bool(text) and text.lower() != historical_default.lower()


def resolve_control_group(
    setting: str | None = "auto",
    moi: str | None = None,
    *,
    perturbo_pool: str | None = None,
    sceptre_group: str | None = None,
) -> dict:
    """Resolve the shared setting into one value per method.

    ``perturbo_pool`` and ``sceptre_group`` are the legacy per-method
    parameters. Either one still on its historical default counts as unset; set
    to anything else it overrides the shared setting for that method alone, and
    the returned log line says so.
    """
    resolved_setting = normalize_setting(setting)
    resolved_moi = normalize_moi(moi)
    notes: list[str] = []

    if resolved_setting == "auto":
        sceptre = SCEPTRE_BY_MOI.get(resolved_moi, "complement")
        # 'from-moi' hands the mapping to the adapter, which also gets to
        # override it from the assignments it can see. Naming a pool here would
        # switch that off.
        perturbo = PERTURBO_HISTORICAL_DEFAULT
        perturbo_effective = PERTURBO_BY_MOI.get(resolved_moi, "auto")
        if resolved_moi == "unknown":
            reason = (
                "auto with no declared low/high MOI: PerTurbo measures the design, "
                "SCEPTRE falls back to 'complement', the only contrast valid for "
                "every design"
            )
            notes.append(reason)
        else:
            reason = f"auto from declared MOI '{resolved_moi}'"
    else:
        sceptre = SCEPTRE_BY_SETTING[resolved_setting]
        perturbo = PERTURBO_BY_SETTING[resolved_setting]
        perturbo_effective = perturbo
        reason = f"{SETTING_PARAM}='{resolved_setting}' set explicitly"

    sceptre_provenance = "auto" if resolved_setting == "auto" else "explicit"
    perturbo_provenance = sceptre_provenance
    overrides: dict[str, str] = {}

    if _is_override(sceptre_group, SCEPTRE_HISTORICAL_DEFAULT):
        sceptre = str(sceptre_group).strip()
        sceptre_provenance = "per-method-override"
        overrides[SCEPTRE_PARAM] = sceptre
    if _is_override(perturbo_pool, PERTURBO_HISTORICAL_DEFAULT):
        perturbo = str(perturbo_pool).strip()
        perturbo_effective = perturbo
        perturbo_provenance = "per-method-override"
        overrides[PERTURBO_PARAM] = perturbo

    # The one combination no amount of configuration can deliver: SCEPTRE has no
    # non-targeting-cell contrast at high MOI. Fail here rather than let the R
    # driver quietly substitute the complement, which is how the two methods
    # ended up answering different questions.
    if sceptre == "nt_cells" and resolved_moi == "high":
        raise ControlGroupError(
            f"control group 'nt_cells' is not available for a high-MOI screen: SCEPTRE's "
            f"non-targeting-cell contrast needs one perturbation per cell. Declared "
            f"Multiplicity_of_infection='high'. Use {SETTING_PARAM}='complement', or "
            f"declare the screen low-MOI if that is what it is."
        )

    perturbo_display = perturbo
    if perturbo != perturbo_effective:
        perturbo_display = f"{perturbo} -> {perturbo_effective}"

    line = (
        f"Control group: declared MOI '{resolved_moi}', {SETTING_PARAM}='{resolved_setting}'"
        f" -> SCEPTRE control_group='{sceptre}', PerTurbo --crt-pool='{perturbo_display}'"
        f" ({reason})."
    )
    if overrides:
        which = ", ".join(f"{name}='{value}'" for name, value in sorted(overrides.items()))
        line += (
            f" {which} overrides the shared setting for that method only:"
            " SCEPTRE and PerTurbo are now deliberately inconsistent and their calls"
            " are not comparable."
        )
    elif sceptre_provenance != "per-method-override" and perturbo_provenance != "per-method-override":
        line += " Both methods contrast against the same cells."

    return {
        "setting": resolved_setting,
        "moi": resolved_moi,
        "sceptre_control_group": sceptre,
        "perturbo_crt_pool": perturbo,
        "perturbo_pool_effective": perturbo_effective,
        "sceptre_provenance": sceptre_provenance,
        "perturbo_provenance": perturbo_provenance,
        "overrides": overrides,
        "reason": reason,
        "notes": notes,
        "log_line": line,
    }


def resolve_from_observations(
    setting: str | None,
    declared_moi: str | None,
    *,
    observed_singleton: bool,
    n_control_only_cells: int,
    n_multi_guide_cells: int = 0,
    perturbo_pool: str | None = None,
    sceptre_group: str | None = None,
) -> dict:
    """The resolution both methods actually receive, settled from the data once.

    ``resolve_control_group`` works from the parameters alone and, for ``auto``,
    hands PerTurbo ``from-moi`` so the adapter can override the declared MOI
    from the assignments it sees. SCEPTRE's driver applied the same rule to its
    own view of the matrix, which is the cis subset, and never changed its
    control group when it did: a screen declared high-MOI whose cells each
    carry one guide ran PerTurbo control-anchored and SCEPTRE against the
    complement. This function makes that decision once, from the full guide
    matrix, and returns explicit values for both methods plus the MOI mode the
    SCEPTRE driver must use, so neither method decides anything on its own.

    ``analysis_moi`` is ``low`` when the screen was declared low or every cell
    carries at most one guide, the rule both drivers already applied; otherwise
    ``high``. Under ``auto`` the non-targeting contrast is chosen only when the
    analysis is low-MOI and at least one cell carries nothing but control guides.
    Explicit settings and per-method overrides pass through as before, except
    that ``from-moi`` is replaced by the pool it maps to, so the adapter's own
    override can no longer fire.
    """
    base = resolve_control_group(
        setting, declared_moi, perturbo_pool=perturbo_pool, sceptre_group=sceptre_group
    )
    declared = base["moi"]
    analysis_moi = "low" if (declared == "low" or observed_singleton) else "high"
    controls_available = int(n_control_only_cells) > 0
    nt_contrast = analysis_moi == "low" and controls_available
    sceptre_from_data = "nt_cells" if nt_contrast else "complement"
    perturbo_from_data = "control-anchored" if nt_contrast else "all-cells"

    sceptre = base["sceptre_control_group"]
    perturbo = base["perturbo_crt_pool"]
    sceptre_provenance = base["sceptre_provenance"]
    perturbo_provenance = base["perturbo_provenance"]
    notes = list(base["notes"])
    if base["setting"] == "auto":
        if sceptre_provenance != "per-method-override":
            sceptre = sceptre_from_data
            sceptre_provenance = "observed"
        if perturbo_provenance != "per-method-override":
            perturbo = perturbo_from_data
            perturbo_provenance = "observed"
        if analysis_moi == "low" and not controls_available:
            notes.append(
                "the analysis is low-MOI but no cell carries only control guides, so "
                "there is no non-targeting pool: both methods use the complement"
            )
        reason = (
            f"auto, settled from the assignments: analysis MOI '{analysis_moi}' "
            f"(declared '{declared}', every cell carries at most one guide: "
            f"{'yes' if observed_singleton else 'no'}), "
            f"{int(n_control_only_cells):,} control-only cell(s)"
        )
    else:
        reason = base["reason"]
    if perturbo == PERTURBO_HISTORICAL_DEFAULT:
        perturbo = perturbo_from_data

    if sceptre == "nt_cells" and analysis_moi != "low":
        raise ControlGroupError(
            "control group 'nt_cells' is not available for this screen: SCEPTRE's "
            "non-targeting-cell contrast needs one perturbation per cell, and the "
            f"assignments carry more than one guide in {int(n_multi_guide_cells):,} cell(s) "
            f"with Multiplicity_of_infection='{declared}'. Use {SETTING_PARAM}='complement'."
        )
    if sceptre == "nt_cells" and not controls_available:
        raise ControlGroupError(
            "control group 'nt_cells' is not available for this screen: no cell carries "
            f"only non-targeting guides, so there is no control population. Use "
            f"{SETTING_PARAM}='complement' or 'auto'."
        )

    overrides = dict(base["overrides"])
    if analysis_moi == "low":
        handling = (
            f"{int(n_multi_guide_cells):,} cell(s) carry more than one guide; SCEPTRE's "
            "low-MOI QC sets them aside and "
            + (
                "PerTurbo's control-anchored pool sets them aside too."
                if perturbo == "control-anchored"
                else f"PerTurbo's '{perturbo}' pool keeps them."
            )
        )
    else:
        handling = "high-MOI analysis: multi-guide cells are kept by both methods."

    line = (
        f"Inference mode: declared MOI '{declared}', observed singleton assignments: "
        f"{'yes' if observed_singleton else 'no'}, control-only cells: {int(n_control_only_cells):,}"
        f" -> analysis MOI '{analysis_moi}', SCEPTRE control_group='{sceptre}', "
        f"PerTurbo --crt-pool='{perturbo}' ({reason})."
    )
    if overrides:
        which = ", ".join(f"{name}='{value}'" for name, value in sorted(overrides.items()))
        line += (
            f" {which} overrides the shared setting for that method only:"
            " SCEPTRE and PerTurbo are now deliberately inconsistent and their calls"
            " are not comparable."
        )
    else:
        line += " Both methods contrast against the same cells."
    line += " " + handling

    return {
        "setting": base["setting"],
        "declared_moi": declared,
        "analysis_moi": analysis_moi,
        "observed_singleton": bool(observed_singleton),
        "n_control_only_cells": int(n_control_only_cells),
        "n_multi_guide_cells": int(n_multi_guide_cells),
        "sceptre_control_group": sceptre,
        "perturbo_crt_pool": perturbo,
        "perturbo_pool_effective": perturbo,
        "sceptre_provenance": sceptre_provenance,
        "perturbo_provenance": perturbo_provenance,
        "overrides": overrides,
        "reason": reason,
        "notes": notes,
        "multi_guide_cell_handling": handling,
        "log_line": line,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--setting", default="auto", help=f"{SETTING_PARAM} value")
    parser.add_argument("--moi", default=None, help="Multiplicity_of_infection value")
    parser.add_argument("--perturbo-pool", default=None, help=f"legacy {PERTURBO_PARAM} value")
    parser.add_argument("--sceptre-group", default=None, help=f"legacy {SCEPTRE_PARAM} value")
    parser.add_argument(
        "--emit",
        default="log",
        choices=["log", "json", "sceptre", "perturbo"],
        help="log line, the whole resolution as JSON, or one method's value",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    resolution = resolve_control_group(
        args.setting,
        args.moi,
        perturbo_pool=args.perturbo_pool,
        sceptre_group=args.sceptre_group,
    )
    if args.emit == "json":
        print(json.dumps(resolution, indent=2))
    elif args.emit == "sceptre":
        print(resolution["sceptre_control_group"])
    elif args.emit == "perturbo":
        print(resolution["perturbo_crt_pool"])
    else:
        print(resolution["log_line"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
