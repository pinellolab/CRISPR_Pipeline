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


def test_final_update_keeps_advanced_execution_dashboard(tmp_path):
    dashboard_dir = tmp_path / "pipeline_dashboard"
    figures = dashboard_dir / "figures"
    svg = dashboard_dir / "svg"
    figures.mkdir(parents=True)
    svg.mkdir()
    # A minimal valid PNG is not required by the renderer; it embeds bytes.
    (figures / "qc.png").write_bytes(b"png-bytes")
    (svg / "plot.svg").write_text("<svg></svg>", encoding="utf-8")
    dashboard = dashboard_dir / "dashboard.html"
    dashboard.write_text(
        '<html><body><img src="svg/plot.svg">'
        '<button data-imgsrc="figures/qc.png">QC</button></body></html>',
        encoding="utf-8",
    )

    trace = tmp_path / "trace.tsv"
    trace.write_text("task_id\tprocess\tstatus\n1\tinputCheck\tCOMPLETED\n", encoding="utf-8")
    output = tmp_path / "telemetry" / "pipeline_execution.html"
    args = type("Args", (), {
        "outdir": tmp_path,
        "run_name": "run-a",
        "trace": trace,
        "source_run_id": "source-a",
        "dashboard_html": output,
        "nextflow_log": tmp_path / "nextflow.log",
        "tail_lines": 30,
        "max_image_bytes": 1_000_000,
    })()

    monitor.render_snapshot(args, "completed", True)
    result = output.read_text(encoding="utf-8")

    assert "Live dependency view" in result
    assert "Pipeline execution" in result
    assert "Pipeline family" in result
    assert "data:image/png;base64," in result
    assert result != dashboard.read_text(encoding="utf-8")


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
