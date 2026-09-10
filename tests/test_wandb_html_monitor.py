import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path


BIN = Path(__file__).parents[1] / "bin"
sys.path.insert(0, str(BIN))
SPEC = importlib.util.spec_from_file_location("wandb_html_monitor", BIN / "wandb_html_monitor.py")
monitor = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(monitor)


def test_discovers_incremental_pipeline_artifacts(tmp_path):
    info = tmp_path / "pipeline_info" / "run-a"
    seqspec = tmp_path / "pipeline_outputs" / "seqspeccheck" / "guide_seqSpec_plots"
    info.mkdir(parents=True)
    seqspec.mkdir(parents=True)
    guide = info / "guide_metadata.validation.json"
    guide.write_text("{}", encoding="utf-8")
    table = seqspec.parent / "guide_position_table.csv"
    table.write_text("Sample,IsWinner\n", encoding="utf-8")
    image = seqspec / "seqSpec_check_plots.png"
    image.write_bytes(b"png")

    paths = monitor.discovered_paths(tmp_path, "run-a")

    assert paths["guide_report"] == guide
    assert paths["seqspec_table"] == table
    assert paths["seqspec_image"] == image


def test_status_file_controls_final_update(tmp_path):
    status = tmp_path / "status.json"
    assert monitor.read_final_status(status) == ("running", False)
    status.write_text(json.dumps({"status": "completed", "exit_code": 0}), encoding="utf-8")
    assert monitor.read_final_status(status) == ("completed", True)


def test_signature_changes_when_trace_changes(tmp_path):
    trace = tmp_path / "trace.tsv"
    status = tmp_path / "status.json"
    trace.write_text("status\n", encoding="utf-8")
    paths = {"missing": tmp_path / "missing"}
    first = monitor.input_signature(paths, trace, status)
    trace.write_text("status\nCOMPLETED\n", encoding="utf-8")
    assert monitor.input_signature(paths, trace, status) != first


def test_wrapper_preserves_nextflow_exit_when_telemetry_is_unavailable(tmp_path):
    fake_nextflow = tmp_path / "nextflow"
    fake_nextflow.write_text("#!/usr/bin/env bash\nexit 7\n", encoding="utf-8")
    fake_nextflow.chmod(0o755)
    state = tmp_path / "state"
    env = os.environ.copy()
    env.pop("WB_IGVF", None)
    env.update({"WANDB_OUTDIR": str(tmp_path / "results"), "WANDB_STATE_DIR": str(state)})

    completed = subprocess.run(
        [str(BIN / "run_with_wandb.sh"), str(fake_nextflow), "run", "main.nf"],
        env=env, text=True, capture_output=True, check=False,
    )

    assert completed.returncode == 7
    assert json.loads((state / "status.json").read_text(encoding="utf-8")) == {
        "status": "failed", "exit_code": 7,
    }
    assert "pipeline continues without W&B" in completed.stderr
