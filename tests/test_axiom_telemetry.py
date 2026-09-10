import importlib.util
import json
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "bin" / "axiom_telemetry.py"
SPEC = importlib.util.spec_from_file_location("axiom_telemetry", MODULE_PATH)
telemetry = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(telemetry)


def test_trace_event_classifies_tools_and_converts_resources():
    event = telemetry.trace_event(
        {"process": "NFCORE_CRISPR:inference_pipeline:inference_perturbo", "name": "fit (1)",
         "status": "COMPLETED", "duration": "1h 2m 3s", "peak_rss": "2.5 GB",
         "%cpu": "125%", "exit": "0", "attempt": "1"},
        {"run_id": "abc"},
    )
    assert event["stage"] == "inference"
    assert event["tool"] == "PerTurbo"
    assert event["duration_ms"] == 3_723_000
    assert event["peak_rss_bytes"] == 2_500_000_000
    assert event["cpu_percent"] == 125


def test_nextflow_timestamp_becomes_rfc3339_utc():
    assert telemetry.event_time("2026-09-10 15:48:09.686") == "2026-09-10T15:48:09.686Z"


def test_dashboard_is_run_scoped_and_layout_matches_charts():
    dashboard = telemetry.dashboard_document("crispr-pipeline", "run-123", "chr8")
    assert dashboard["refreshTime"] == 60
    assert dashboard["owner"] == "X-AXIOM-EVERYONE"
    assert {item["i"] for item in dashboard["layout"]} == {item["id"] for item in dashboard["charts"]}
    queries = [chart["query"]["apl"] for chart in dashboard["charts"] if "query" in chart]
    assert all('run_id == "run-123"' in query for query in queries)
    assert {"Input and guide QC", "SeqSpec QC by sample", "QC images and reports generated"} <= {
        chart["name"] for chart in dashboard["charts"]
    }


def test_trace_reader_deduplicates_terminal_rows(tmp_path):
    trace = tmp_path / "trace.tsv"
    trace.write_text("task_id\thash\tattempt\tprocess\tstatus\n1\taa\t1\tp\tCOMPLETED\n", encoding="utf-8")
    seen = set()
    assert len(telemetry.read_new_trace_rows(trace, seen)) == 1
    assert telemetry.read_new_trace_rows(trace, seen) == []


def test_hook_reader_preserves_run_scope_and_offset(tmp_path):
    hooks = tmp_path / "hooks.jsonl"
    hooks.write_text(json.dumps({"event_type": "workflow_hook_complete", "status": "SUCCEEDED"}) + "\n")
    events, offset = telemetry.read_new_hook_events(hooks, 0, {"run_id": "expected"})
    assert events[0]["run_id"] == "expected"
    assert telemetry.read_new_hook_events(hooks, offset, {"run_id": "expected"})[0] == []


def test_qc_metric_flattening_is_numeric_and_bounded(tmp_path):
    source = tmp_path / "pipeline_qc_metrics.json"
    source.write_text(json.dumps({"cells": 10, "nested": {"rate": 0.5}, "label": "skip"}), encoding="utf-8")
    events = telemetry.flatten_qc_metrics(source, {"run_id": "abc"})
    assert {(event["metric"], event["value"]) for event in events} == {("cells", 10), ("nested.rate", 0.5)}


def test_dependency_events_extract_graphviz_edges(tmp_path):
    dag = tmp_path / "dag.dot"
    dag.write_text('digraph x {\np0 [label="mappingGuide"];\np1 [label="CreateMuData"];\np0 -> p1;\n}\n')
    events = telemetry.dependency_events(dag, {"run_id": "abc"})
    assert [(event["upstream_process"], event["downstream_process"]) for event in events] == [
        ("mappingGuide", "CreateMuData")
    ]


def test_live_qc_discovers_inputs_seqspec_and_image_inventory(tmp_path):
    run_name = "chr8"
    info = tmp_path / "pipeline_info" / run_name
    seqspec = tmp_path / "pipeline_outputs" / "seqspeccheck"
    plots = seqspec / "guide_seqSpec_plots"
    info.mkdir(parents=True)
    plots.mkdir(parents=True)
    (tmp_path / "pipeline_info" / "original_samplesheet.csv").write_text(
        "file_modality,measurement_sets\nscRNA,set1\ngRNA,set1\nscRNA,set2\n", encoding="utf-8"
    )
    (info / "guide_metadata.validation.json").write_text(
        json.dumps({"row_count": 12, "control_rows": 2, "valid": True}), encoding="utf-8"
    )
    (seqspec / "guide_position_table.csv").write_text(
        "Sample,Config,TotalHits,HitRatio,PosPurity,FlankPurity,Gini,FinalScore,IsWinner\n"
        "s1,R2_Fwd,100,0.1,0.9,0.8,0.7,3.4,True\n", encoding="utf-8"
    )
    (plots / "seqSpec_check_plots.png").write_bytes(b"png")
    seen = {}
    events = telemetry.discover_live_qc(tmp_path, run_name, {"run_id": "abc"}, seen)
    metrics = {(event.get("metric"), event.get("value")) for event in events if event["event_type"] == "qc_metric"}
    assert ("input.samplesheet_rows", 3) in metrics
    assert ("guide_metadata.row_count", 12) in metrics
    assert ("seqspec.HitRatio", 0.1) in metrics
    assert any(event.get("artifact_name") == "seqSpec_check_plots.png" for event in events)
    assert telemetry.discover_live_qc(tmp_path, run_name, {"run_id": "abc"}, seen) == []


def test_lifecycle_property_access_is_inside_fail_open_guard():
    source = (Path(__file__).parents[1] / "subworkflows/local/utils_nfcore_crispr_pipeline/main.nf").read_text()
    complete = source[source.index("workflow.onComplete"):source.index("workflow.onError")]
    error = source[source.index("workflow.onError"):source.index("// The sidecar tails")]
    assert complete.index("try {") < complete.index("workflow?.success") < complete.index("catch (Exception error)")
    assert error.index("try {") < error.index("workflow?.errorMessage") < error.index("catch (Exception error)")
