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
