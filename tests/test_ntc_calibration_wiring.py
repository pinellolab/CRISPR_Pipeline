"""The non-targeting calibration pass is wired the way the scripts expect.

Static checks on the Nextflow sources, the configuration and the schema, plus
a parse of every flag each module hands to its script against that script's
own parser, so a renamed option cannot slip through to a run.
"""

import json
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

PARAMS = (
    "ENABLE_NTC_CALIBRATION",
    "NTC_CALIBRATION_max_local_pairs",
    "NTC_CALIBRATION_max_global_pairs",
    "NTC_CALIBRATION_seed",
)
SCEPTRE_PROCESSES = ("sceptre_chunk_prepare", "inference_sceptre", "sceptre_chunk_merge")


def _read(relative):
    return (REPO_ROOT / relative).read_text()


def _flags(script_text):
    return set(re.findall(r"(--[a-z][a-z0-9-]*)", script_text))


def _parser_flags(module_name):
    module = __import__(module_name)
    return {opt for action in module.build_parser()._actions for opt in action.option_strings}


def test_the_pass_is_called_from_both_two_method_branches_and_emitted():
    workflow = _read("subworkflows/local/inference_pipeline/main.nf")
    assert "include { ntc_calibration_pipeline } from '../ntc_calibration_pipeline'" in workflow
    assert workflow.count("ntc_calibration_pipeline(") == 2
    assert "calibration_output = calibration_output" in workflow
    # The flag is decided on its text: Nextflow 26 passes `--X false` as 'false'.
    assert "calibration_flag.toString().trim().toLowerCase() in ['true', '1', 'yes']" in workflow
    assert "assets/calibration_empty" in workflow
    assert (REPO_ROOT / "assets/calibration_empty/.gitkeep").exists()


def test_the_subworkflow_reuses_the_sceptre_processes_under_their_own_names():
    subworkflow = _read("subworkflows/local/ntc_calibration_pipeline/main.nf")
    for name in SCEPTRE_PROCESSES:
        assert f"include {{ {name} }} from" in subworkflow, name
        assert f"{name} as " not in subworkflow, f"{name} must not be aliased"
    assert "JsonSlurper().parseText(f.text).status" in subworkflow
    assert "ntc_calibration_skipped(" in subworkflow
    assert "calibration_output = Report.calibration_output.mix(Skipped.calibration_output)" in subworkflow


def test_both_workflow_branches_fold_the_calibration_into_the_benchmark_outputs():
    workflow = _read("workflows/crispr_pipeline/main.nf")
    assert "include { collect_benchmark_outputs }" in workflow
    assert workflow.count("collect_benchmark_outputs(benchmark_output_dir, Inference.calibration_output).collected") == 2


def test_parameters_are_declared_in_the_config_and_the_schema():
    config = _read("nextflow.config")
    schema = json.loads(_read("nextflow_schema.json"))
    declared = {}
    for group in schema.get("definitions", schema.get("$defs", {})).values():
        declared.update(group.get("properties", {}))
    for name in PARAMS:
        assert re.search(rf"^\s*{name} = ", config, re.M), name
        assert name in declared, name
    assert declared["ENABLE_NTC_CALIBRATION"]["default"] is True


def test_the_new_processes_have_resource_entries():
    config = _read("nextflow.config")
    for name in ("resolve_inference_mode", "prepare_ntc_calibration", "ntc_calibration_report", "ntc_calibration_skipped", "collect_benchmark_outputs"):
        assert re.search(rf"withName: '[^']*\b{name}\b[^']*'", config), name


@pytest.mark.parametrize(
    ("module", "script"),
    [
        ("modules/local/prepare_ntc_calibration/main.nf", "prepare_ntc_calibration"),
        ("modules/local/ntc_calibration_report/main.nf", "ntc_calibration_report"),
        ("modules/local/resolve_inference_mode/main.nf", "resolve_inference_mode"),
    ],
)
def test_every_flag_a_module_passes_exists_on_its_script(module, script):
    passed = _flags(_read(module))
    known = _parser_flags(script)
    unknown = sorted(passed - known)
    assert not unknown, f"{module} passes flags {script} does not define: {unknown}"


def test_nextflow_resolves_the_parameters():
    nextflow = shutil.which("nextflow") or str(pathlib.Path.home() / ".local/bin/nextflow")
    if not pathlib.Path(nextflow).exists():
        pytest.skip("nextflow is not installed")
    result = subprocess.run([nextflow, "config", str(REPO_ROOT), "-flat"], capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    for name in PARAMS:
        assert f"params.{name} = " in result.stdout, name
