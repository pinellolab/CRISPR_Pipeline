"""The one setting that decides which cells a perturbation is compared against.

The resolver is the only place the SCEPTRE/PerTurbo vocabulary is translated, so
these tests pin the whole table, the precedence of the older per-method
parameters, and the two configurations that must fail rather than be silently
substituted. The Groovy mirror the pipeline actually runs is checked against the
Python table at the bottom, and the module argv is checked to confirm both
processes are driven by the resolution rather than by the raw params.
"""

import pathlib
import re
import shutil
import subprocess
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import control_group as cg


# --- the resolution table: every setting against every MOI ------------------


@pytest.mark.parametrize(
    ("setting", "moi", "sceptre", "perturbo", "perturbo_effective"),
    [
        # 'auto' takes the declared MOI. These two rows are the pre-existing
        # behaviour for PerTurbo, so a run on the historical defaults is unchanged.
        ("auto", "low", "nt_cells", "from-moi", "control-anchored"),
        ("auto", "high", "complement", "from-moi", "all-cells"),
        # Neither: PerTurbo measures the design, SCEPTRE cannot, so it takes the
        # contrast that is valid for every design.
        ("auto", None, "complement", "from-moi", "auto"),
        ("auto", "", "complement", "from-moi", "auto"),
        ("auto", "medium", "complement", "from-moi", "auto"),
        # Named outright: the MOI does not enter into it.
        ("nt_cells", "low", "nt_cells", "control-anchored", "control-anchored"),
        ("nt_cells", None, "nt_cells", "control-anchored", "control-anchored"),
        ("complement", "low", "complement", "all-cells", "all-cells"),
        ("complement", "high", "complement", "all-cells", "all-cells"),
        ("complement", None, "complement", "all-cells", "all-cells"),
    ],
)
def test_resolution_table(setting, moi, sceptre, perturbo, perturbo_effective):
    resolution = cg.resolve_control_group(setting, moi)

    assert resolution["sceptre_control_group"] == sceptre
    assert resolution["perturbo_crt_pool"] == perturbo
    assert resolution["perturbo_pool_effective"] == perturbo_effective


def test_auto_keeps_perturbos_own_mapping_so_the_adapter_can_override_it():
    """'auto' must pass 'from-moi', not the mapped pool.

    The adapter overrides a declared MOI when the assignments disagree with it --
    one perturbation per cell plus control-only cells means the control-anchored
    pool whatever the samplesheet says -- and it only does that when the pool was
    handed to it as 'from-moi'. Naming the pool here would switch that off.
    """
    for moi in ("low", "high", None):
        assert cg.resolve_control_group("auto", moi)["perturbo_crt_pool"] == "from-moi"


@pytest.mark.parametrize("setting", ["AUTO", " nt_cells ", "Complement"])
def test_setting_is_case_and_whitespace_insensitive(setting):
    assert cg.resolve_control_group(setting, "low")["setting"] == setting.strip().lower()


def test_unset_setting_is_auto():
    assert cg.resolve_control_group(None, "low")["sceptre_control_group"] == "nt_cells"


# --- rejected values --------------------------------------------------------


@pytest.mark.parametrize("setting", ["nt", "control-anchored", "all-cells", "from-moi", "none", ""])
def test_invalid_setting_names_the_accepted_values(setting):
    with pytest.raises(cg.ControlGroupError) as excinfo:
        cg.resolve_control_group(setting, "low")

    message = str(excinfo.value)
    assert "INFERENCE_control_group" in message
    for accepted in ("auto", "nt_cells", "complement"):
        assert accepted in message


# --- the one combination SCEPTRE cannot provide -----------------------------


def test_nt_cells_at_high_moi_is_refused_not_substituted():
    with pytest.raises(cg.ControlGroupError) as excinfo:
        cg.resolve_control_group("nt_cells", "high")

    message = str(excinfo.value)
    assert "nt_cells" in message
    assert "high-MOI" in message
    assert "complement" in message


def test_nt_cells_at_high_moi_is_refused_through_the_per_method_override_too():
    with pytest.raises(cg.ControlGroupError):
        cg.resolve_control_group("auto", "high", sceptre_group="nt_cells")


def test_high_moi_auto_is_fine():
    assert cg.resolve_control_group("auto", "high")["sceptre_control_group"] == "complement"


# --- precedence of the older per-method parameters --------------------------


def test_historical_defaults_count_as_unset():
    """Every config written before the shared setting existed carries these."""
    shared = cg.resolve_control_group("nt_cells", "low")
    legacy = cg.resolve_control_group(
        "nt_cells", "low", perturbo_pool="from-moi", sceptre_group="complement"
    )

    assert legacy["sceptre_control_group"] == shared["sceptre_control_group"] == "nt_cells"
    assert legacy["perturbo_crt_pool"] == shared["perturbo_crt_pool"] == "control-anchored"
    assert legacy["sceptre_provenance"] == "explicit"
    assert legacy["perturbo_provenance"] == "explicit"
    assert legacy["overrides"] == {}


def test_perturbo_override_wins_for_perturbo_only():
    resolution = cg.resolve_control_group("nt_cells", "low", perturbo_pool="all-cells")

    assert resolution["perturbo_crt_pool"] == "all-cells"
    assert resolution["perturbo_provenance"] == "per-method-override"
    # SCEPTRE is untouched.
    assert resolution["sceptre_control_group"] == "nt_cells"
    assert resolution["sceptre_provenance"] == "explicit"
    assert resolution["overrides"] == {"INFERENCE_PERTURBO_CRT_POOL": "all-cells"}


def test_sceptre_override_wins_for_sceptre_only():
    resolution = cg.resolve_control_group("complement", "low", sceptre_group="nt_cells")

    assert resolution["sceptre_control_group"] == "nt_cells"
    assert resolution["sceptre_provenance"] == "per-method-override"
    assert resolution["perturbo_crt_pool"] == "all-cells"
    assert resolution["perturbo_provenance"] == "explicit"
    assert resolution["overrides"] == {"INFERENCE_SCEPTRE_control_group": "nt_cells"}


def test_both_overrides_are_reported_together():
    resolution = cg.resolve_control_group(
        "auto", "low", perturbo_pool="all-cells", sceptre_group="nt_cells"
    )

    assert resolution["overrides"] == {
        "INFERENCE_PERTURBO_CRT_POOL": "all-cells",
        "INFERENCE_SCEPTRE_control_group": "nt_cells",
    }


def test_an_override_says_the_methods_are_inconsistent():
    line = cg.resolve_control_group("nt_cells", "low", perturbo_pool="all-cells")["log_line"]

    assert "INFERENCE_PERTURBO_CRT_POOL='all-cells'" in line
    assert "overrides the shared setting for that method only" in line
    assert "deliberately inconsistent" in line
    assert "Both methods contrast against the same cells" not in line


# --- the log line -----------------------------------------------------------


def test_log_line_names_the_moi_the_setting_both_methods_and_the_reason():
    line = cg.resolve_control_group("auto", "low")["log_line"]

    assert line.count("\n") == 0
    assert "declared MOI 'low'" in line
    assert "INFERENCE_control_group='auto'" in line
    assert "SCEPTRE control_group='nt_cells'" in line
    assert "PerTurbo --crt-pool='from-moi -> control-anchored'" in line
    assert "auto from declared MOI 'low'" in line
    assert "Both methods contrast against the same cells." in line


def test_log_line_explains_the_undeclared_moi_fallback():
    line = cg.resolve_control_group("auto", None)["log_line"]

    assert "declared MOI 'unknown'" in line
    assert "PerTurbo measures the design" in line
    assert "valid for every design" in line


# --- the Groovy mirror the pipeline actually runs ---------------------------


GROOVY = (REPO_ROOT / "modules/local/control_group/main.nf").read_text()


def _groovy_map(function_name):
    body = re.search(rf"def {function_name}\(\) \{{ return \[(.*?)\] }}", GROOVY).group(1)
    return dict(re.findall(r"'([^']+)':\s*'([^']+)'", body))


def test_groovy_mirror_agrees_with_the_python_table():
    assert _groovy_map("controlGroupSceptreBySetting") == cg.SCEPTRE_BY_SETTING
    assert _groovy_map("controlGroupPerturboBySetting") == cg.PERTURBO_BY_SETTING
    assert _groovy_map("controlGroupPerturboByMoi") == cg.PERTURBO_BY_MOI
    assert _groovy_map("controlGroupSceptreByMoi") == cg.SCEPTRE_BY_MOI


def test_groovy_mirror_agrees_on_the_accepted_values_and_historical_defaults():
    settings = re.search(r"def controlGroupSettings\(\) \{\s*return \[(.*?)\]", GROOVY, re.S).group(1)
    assert tuple(re.findall(r"'([^']+)'", settings)) == cg.SETTINGS
    assert f"controlGroupPerturboHistoricalDefault() {{ return '{cg.PERTURBO_HISTORICAL_DEFAULT}' }}" in GROOVY
    assert f"controlGroupSceptreHistoricalDefault() {{ return '{cg.SCEPTRE_HISTORICAL_DEFAULT}' }}" in GROOVY


def test_groovy_mirror_refuses_nt_cells_at_high_moi():
    assert "sceptre == 'nt_cells' && resolvedMoi == 'high'" in GROOVY
    assert "is not available for a high-MOI screen" in GROOVY


# --- the pipeline wiring ----------------------------------------------------


def test_the_shared_setting_is_declared_with_auto_as_the_default():
    config = (REPO_ROOT / "nextflow.config").read_text()

    assert re.search(r"^\s*INFERENCE_control_group = 'auto'$", config, re.M)
    # The per-method parameters stay on their historical values, which the
    # resolver reads as "not set".
    assert re.search(r"^\s*INFERENCE_PERTURBO_CRT_POOL = 'from-moi'$", config, re.M)
    assert re.search(r"^\s*INFERENCE_SCEPTRE_control_group = 'complement'$", config, re.M)


def test_the_resolution_happens_once_in_the_inference_subworkflow():
    workflow = (REPO_ROOT / "subworkflows/local/inference_pipeline/main.nf").read_text()

    assert "include { resolveControlGroupFromParams }" in workflow
    assert workflow.count("resolveControlGroupFromParams()") == 1
    assert "log.info(control_group.log_line)" in workflow


def test_both_processes_are_driven_by_the_resolution_not_by_the_raw_params():
    sceptre_module = (REPO_ROOT / "modules/local/inference_sceptre/main.nf").read_text()
    perturbo_module = (REPO_ROOT / "modules/local/inference_perturbo/main.nf").read_text()
    workflow = (REPO_ROOT / "subworkflows/local/inference_pipeline/main.nf").read_text()

    # SCEPTRE's control group reaches the R driver from the resolution.
    assert "${params.INFERENCE_SCEPTRE_control_group}" not in sceptre_module
    assert "${control_group}" in sceptre_module
    # PerTurbo's pool likewise, with the provenance for the recorded metadata.
    assert "${params.INFERENCE_PERTURBO_CRT_POOL}" not in perturbo_module
    assert "--crt-pool ${control_group.perturbo_crt_pool}" in perturbo_module
    assert "--control-group-setting ${control_group.setting}" in perturbo_module
    assert "--control-group-provenance ${control_group.perturbo_provenance}" in perturbo_module

    # Every call site passes it: three SCEPTRE calls, three PerTurbo calls.
    assert workflow.count("control_group.sceptre_control_group)") == 3
    assert len(re.findall(r"inference_perturbo\([^)]*control_group\)", workflow)) == 3


def test_the_adapter_accepts_the_provenance_arguments():
    parser = __import__("perturbo_v2_pipeline_adapter").build_parser()
    args = parser.parse_args(
        [
            "--input", "in.h5mu",
            "--per-element-output", "e.tsv.gz",
            "--per-guide-output", "g.tsv.gz",
            "--crt-pool", "control-anchored",
            "--control-group-setting", "nt_cells",
            "--control-group-provenance", "explicit",
        ]
    )

    assert args.crt_pool == "control-anchored"
    assert args.control_group_setting == "nt_cells"
    assert args.control_group_provenance == "explicit"


# --- the SCEPTRE driver ----------------------------------------------------


SCEPTRE_R = (REPO_ROOT / "bin/inference_sceptre.R").read_text()


def test_the_r_driver_no_longer_forces_the_complement():
    """The old driver assigned 'complement' unconditionally and warned.

    Behaviour is now: honour what was asked, decide from the object when nothing
    was asked, and stop when what was asked cannot be delivered.
    """
    assert 'Overriding control_group=' not in SCEPTRE_R
    assert 'args_list$control_group <- "complement"' not in SCEPTRE_R
    assert "args_list$control_group <- control_group_resolution$control_group" in SCEPTRE_R


def test_the_r_driver_refuses_nt_cells_at_high_moi():
    assert 'requested_control_group == "nt_cells" && !low_moi' in SCEPTRE_R
    assert "this is a high-MOI analysis" in SCEPTRE_R
    # And refuses it when there is no control population either.
    assert 'requested_control_group == "nt_cells" && !nt_cells_available' in SCEPTRE_R


def test_the_r_driver_treats_auto_as_no_request():
    assert 'tolower(raw) %in% c("auto", "default")' in SCEPTRE_R


def test_the_r_driver_rejects_an_unknown_control_group():
    assert 'accepted_control_groups <- c("nt_cells", "complement")' in SCEPTRE_R
    assert "is not a SCEPTRE control group" in SCEPTRE_R


def test_the_r_resolution_is_a_function_of_its_own_so_it_can_be_checked():
    assert "resolve_sceptre_control_group <- function(requested, is_low_moi, n_nt_cells)" in SCEPTRE_R
    assert "control_group_resolution <- resolve_sceptre_control_group(" in SCEPTRE_R


def test_the_r_driver_passes_its_own_checks():
    """Run the R-level check of the extracted resolver against a real R.

    bin/inference_sceptre.R cannot be sourced without sceptre and the
    Bioconductor stack, so tests/control_group_sceptre_check.R pulls out the two
    self-contained functions and exercises them. Skipped where R is absent.
    """
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")

    check = REPO_ROOT / "tests/control_group_sceptre_check.R"
    result = subprocess.run(
        [rscript, str(check), str(REPO_ROOT / "bin/inference_sceptre.R")],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "all checks passed" in result.stdout


# --- the recorded provenance ----------------------------------------------


def test_the_resolution_is_recorded_beside_the_results():
    adapter = (BIN_DIR / "perturbo_v2_pipeline_adapter.py").read_text()

    assert "control_group_resolution.json" in adapter
    # PerTurbo's own crt_metadata.json names the pool but cannot know where it
    # came from; the adapter adds the pipeline's side to each copy.
    assert 'rglob("crt_metadata.json")' in adapter
    assert '"pipeline_control_group"' in adapter
    # The SCEPTRE side writes the same facts.
    assert "sceptre_control_group.json" in SCEPTRE_R
    assert "control_group_requested" in SCEPTRE_R
    sceptre_module = (REPO_ROOT / "modules/local/inference_sceptre/main.nf").read_text()
    assert 'path "sceptre_control_group.json", optional: true' in sceptre_module


def test_the_adapter_records_pool_provenance_and_reason():
    import perturbo_v2_pipeline_adapter as adapter

    source = (BIN_DIR / "perturbo_v2_pipeline_adapter.py").read_text()
    for key in ('"crt_pool"', '"crt_pool_requested"', '"declared_moi"', '"provenance"', '"reason"'):
        assert key in source
    assert adapter.SETTING_PARAM == "INFERENCE_control_group"


def test_record_control_group_annotates_every_crt_metadata_copy(tmp_path):
    import json

    import perturbo_v2_pipeline_adapter as adapter

    for name in ("element", "guide"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "crt_metadata.json").write_text(json.dumps({"pool": "control-anchored"}))
    (tmp_path / "element" / "not_metadata.json").write_text("{}")

    record = {"method": "perturbo", "crt_pool": "control-anchored", "provenance": "auto"}
    adapter._record_control_group(tmp_path, record)

    assert json.loads((tmp_path / "control_group_resolution.json").read_text()) == record
    for name in ("element", "guide"):
        metadata = json.loads((tmp_path / name / "crt_metadata.json").read_text())
        assert metadata["pool"] == "control-anchored"
        assert metadata["pipeline_control_group"] == record


def test_record_control_group_survives_a_malformed_metadata_file(tmp_path):
    import perturbo_v2_pipeline_adapter as adapter

    (tmp_path / "crt_metadata.json").write_text("not json")
    adapter._record_control_group(tmp_path, {"crt_pool": "all-cells"})

    assert (tmp_path / "control_group_resolution.json").exists()
    assert (tmp_path / "crt_metadata.json").read_text() == "not json"


def test_record_control_group_does_nothing_without_a_record(tmp_path):
    import perturbo_v2_pipeline_adapter as adapter

    adapter._record_control_group(tmp_path, None)

    assert list(tmp_path.iterdir()) == []


# --- the documented behaviour change --------------------------------------


def test_the_changelog_says_low_moi_sceptre_results_change():
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text()

    assert "INFERENCE_control_group" in changelog
    # Pin the claim, not one phrasing of it: the changelog must say somewhere
    # that low-MOI SCEPTRE output is not comparable to earlier runs. Matching a
    # literal sentence made this fail on an editorial rewrite that kept the warning.
    warning = re.search(
        r"low[- ]MOI[^.]{0,400}?(differs?|will differ|not comparable|change[sd]?)",
        changelog,
        re.IGNORECASE | re.DOTALL,
    )
    assert warning, "CHANGELOG.md must warn that low-MOI SCEPTRE results change"
