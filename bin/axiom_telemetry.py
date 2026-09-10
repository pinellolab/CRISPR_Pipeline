#!/usr/bin/env python3
"""Best-effort, size-bounded Axiom telemetry for a Nextflow execution.

The runner preserves the wrapped Nextflow exit code. Network, dashboard, parsing,
or telemetry failures are warnings only and can never fail the scientific run.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


DEFAULT_MAX_BYTES = 20_000_000
DEFAULT_MAX_EVENT_BYTES = 8_192
TERMINAL_STATES = {"COMPLETED", "FAILED", "CACHED", "ABORTED"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def event_time(value: Any) -> str:
    """Normalize Nextflow's local trace timestamp to RFC3339 UTC when possible."""
    if not value:
        return utc_now()
    text = str(value).strip()
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")
    except ValueError:
        return utc_now()


def safe_text(value: Any, limit: int = 2048) -> str:
    text = "" if value is None else str(value)
    text = re.sub(r"(?i)(authorization:\s*bearer\s+|xaat-)[^\s'\"]+", r"\1<redacted>", text)
    return text[:limit]


def process_metadata(process_name: str) -> dict[str, str]:
    leaf = process_name.split(":")[-1]
    low = leaf.lower()
    rules = [
        (("mappingguide", "mappinghashing", "mappingscrna"), "mapping", "kallisto|bustools"),
        (("seqspec",), "input-qc", "seqspec"),
        (("guide_assignment_sceptre",), "guide-assignment", "SCEPTRE"),
        (("guide_assignment_cleanser",), "guide-assignment", "CRISPRcleanR"),
        (("perturbo",), "inference", "PerTurbo"),
        (("inference_sceptre", "sceptre_chunk"), "inference", "SCEPTRE"),
        (("preprocess", "doublets", "filter"), "preprocessing", "Scanpy"),
        (("dashboard", "additional_qc", "evaluation"), "reporting", "pipeline dashboard"),
        (("catalog", "mergedresults", "mergemudata"), "results", "MuData|Parquet"),
    ]
    for needles, stage, tool in rules:
        if any(needle in low for needle in needles):
            return {"process_leaf": leaf, "stage": stage, "tool": tool}
    return {"process_leaf": leaf, "stage": "pipeline", "tool": leaf}


def trace_event(row: dict[str, str], base: dict[str, Any]) -> dict[str, Any]:
    process_name = row.get("process") or row.get("name") or "unknown"
    event = dict(base)
    event.update(process_metadata(process_name))
    event.update(
        {
            "_time": event_time(row.get("complete") or row.get("start")),
            "event_type": "process_completed",
            "message": f"{process_name} {row.get('status', 'UNKNOWN')}",
            "process": process_name,
            "task_name": row.get("name", ""),
            "task_id": row.get("task_id", ""),
            "task_hash": row.get("hash", ""),
            "native_id": row.get("native_id", ""),
            "status": row.get("status", "UNKNOWN"),
            "exit_code": to_number(row.get("exit")),
            "attempt": to_number(row.get("attempt")),
            "duration_ms": duration_ms(row.get("realtime") or row.get("duration")),
            "cpu_percent": to_number((row.get("%cpu") or "").rstrip("%")),
            "peak_rss_bytes": memory_bytes(row.get("peak_rss") or row.get("rss")),
            "read_bytes": memory_bytes(row.get("read_bytes") or row.get("rchar")),
            "write_bytes": memory_bytes(row.get("write_bytes") or row.get("wchar")),
            "container": safe_text(row.get("container", ""), 512),
            "cached": row.get("status") == "CACHED",
        }
    )
    return event


def to_number(value: Any) -> int | float | None:
    if value in (None, "", "-"):
        return None
    try:
        number = float(str(value).replace(",", ""))
        return int(number) if number.is_integer() else number
    except ValueError:
        return None


def memory_bytes(value: Any) -> int | None:
    if value in (None, "", "-"):
        return None
    match = re.fullmatch(r"\s*([0-9.]+)\s*([KMGTPE]?B)?\s*", str(value), re.I)
    if not match:
        return to_number(value)  # type: ignore[return-value]
    scale = {"": 1, "B": 1, "KB": 1000, "MB": 1000**2, "GB": 1000**3,
             "TB": 1000**4, "PB": 1000**5, "EB": 1000**6}
    return int(float(match.group(1)) * scale[(match.group(2) or "").upper()])


def duration_ms(value: Any) -> int | None:
    if value in (None, "", "-"):
        return None
    text = str(value).strip().lower()
    if text.isdigit():
        return int(text)
    total = 0.0
    matches = re.findall(r"([0-9.]+)\s*(ms|us|s|m|h|d)", text)
    if not matches:
        return None
    factors = {"us": 0.001, "ms": 1, "s": 1000, "m": 60_000, "h": 3_600_000, "d": 86_400_000}
    for amount, unit in matches:
        total += float(amount) * factors[unit]
    return int(total)


class AxiomClient:
    def __init__(self, dataset: str, token_env: str, ingest_url: str, api_url: str,
                 timeout: float, max_bytes: int, max_event_bytes: int, state_dir: Path):
        self.dataset = dataset
        self.token_env = token_env
        self.ingest_url = ingest_url.rstrip("/")
        self.api_url = api_url.rstrip("/")
        self.timeout = timeout
        self.max_bytes = max_bytes
        self.max_event_bytes = max_event_bytes
        self.sent_bytes = 0
        self.dropped_events = 0
        self.disabled_reason = ""
        self.state_dir = state_dir
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.warning_log = self.state_dir / "axiom_telemetry.warnings.log"

    @property
    def token(self) -> str:
        return os.environ.get(self.token_env, "")

    def warn(self, message: str) -> None:
        line = f"{utc_now()} {safe_text(message, 2048)}\n"
        try:
            if self.warning_log.exists() and self.warning_log.stat().st_size >= 1_000_000:
                return
            with self.warning_log.open("a", encoding="utf-8") as handle:
                handle.write(line)
        except OSError:
            pass
        print(f"WARN: Axiom telemetry: {safe_text(message, 512)}", file=sys.stderr)

    def request(self, method: str, url: str, payload: Any | None = None) -> tuple[int, str]:
        if not self.token:
            raise RuntimeError(f"environment variable {self.token_env} is not set")
        body = None if payload is None else json.dumps(payload, separators=(",", ":")).encode()
        request = urllib.request.Request(
            url, data=body, method=method,
            headers={"Authorization": f"Bearer {self.token}", "Content-Type": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return response.status, response.read(262_144).decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            detail = exc.read(8_192).decode("utf-8", "replace")
            raise RuntimeError(f"HTTP {exc.code}: {safe_text(detail)}") from exc

    def ingest(self, events: Iterable[dict[str, Any]]) -> bool:
        if self.disabled_reason:
            return False
        accepted = []
        for event in events:
            encoded = json.dumps(event, separators=(",", ":"), default=str).encode()
            if len(encoded) > self.max_event_bytes:
                event = dict(event)
                event["message"] = safe_text(event.get("message", ""), 512)
                event["telemetry_truncated"] = True
                encoded = json.dumps(event, separators=(",", ":"), default=str).encode()
            if len(encoded) > self.max_event_bytes or self.sent_bytes + len(encoded) > self.max_bytes:
                self.dropped_events += 1
                continue
            accepted.append(event)
            self.sent_bytes += len(encoded)
        if not accepted:
            return False
        dataset = urllib.parse.quote(self.dataset, safe="")
        if "{dataset}" in self.ingest_url:
            url = self.ingest_url.replace("{dataset}", dataset)
        else:
            url = f"{self.ingest_url}/{dataset}"
        try:
            self.request("POST", url, accepted)
            return True
        except Exception as exc:  # telemetry is deliberately fail-open
            self.warn(f"ingest failed; pipeline continues: {exc}")
            return False

    def create_dashboard(self, run_id: str, run_name: str) -> str | None:
        dashboard = dashboard_document(self.dataset, run_id, run_name)
        uid = f"crispr-pipeline-{run_id.lower()}"
        payload = {"dashboard": dashboard, "uid": uid, "overwrite": True,
                   "message": f"Provision dashboard for Nextflow run {run_id}"}
        try:
            self.request("POST", f"{self.api_url}/v2/dashboards", payload)
            return uid
        except Exception as exc:
            self.warn(f"dashboard provisioning failed; pipeline continues: {exc}")
            return None


def dashboard_document(dataset: str, run_id: str, run_name: str) -> dict[str, Any]:
    ds = dataset.replace("'", "")
    rid = run_id.replace('"', "")
    prefix = f"['{ds}'] | where run_id == \"{rid}\""
    charts = [
        {"id": "run-status", "name": "Latest run status", "type": "Table",
         "query": {"apl": prefix + " | where event_type in (\"run_started\", \"run_completed\") | top 1 by _time desc | project _time, status, run_name, elapsed_seconds"}},
        {"id": "completed-tasks", "name": "Completed tasks", "type": "Statistic", "unit": "Abbreviated",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" | summarize count()"}},
        {"id": "failed-tasks", "name": "Failed tasks", "type": "Statistic", "unit": "Abbreviated", "colorScheme": "Red",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" and status == \"FAILED\" | summarize count()"}},
        {"id": "progress", "name": "Task completions over time", "type": "TimeSeries",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" | summarize tasks=count() by bin_auto(_time), status"}},
        {"id": "stage-status", "name": "Process status by pipeline stage", "type": "Table",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" | summarize tasks=count(), runtime_seconds=sum(duration_ms) / 1000 by stage, status | order by stage asc"}},
        {"id": "slowest", "name": "Slowest processes", "type": "Table", "customUnits": "s",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" | summarize tasks=count(), mean_seconds=avg(duration_ms) / 1000, max_seconds=max(duration_ms) / 1000 by process_leaf, tool | top 15 by max_seconds desc"}},
        {"id": "resources", "name": "Peak memory by process", "type": "Table", "customUnits": "GB",
         "query": {"apl": prefix + " | where event_type == \"process_completed\" | summarize peak_gb=max(peak_rss_bytes) / 1000000000, cpu_percent=avg(cpu_percent) by process_leaf | top 15 by peak_gb desc"}},
        {"id": "qc-metrics", "name": "Pipeline QC metrics", "type": "Table",
         "query": {"apl": prefix + " | where event_type == \"qc_metric\" | project metric, value, unit, source | take 200"}},
        {"id": "dependencies", "name": "Process dependencies", "type": "Table",
         "query": {"apl": prefix + " | where event_type == \"process_dependency\" | distinct upstream_process, downstream_process | take 200"}},
        {"id": "event-log", "name": "Run event log", "type": "LogStream",
         "query": {"apl": prefix + " | project _time, event_type, status, stage, process_leaf, message | take 300"}},
    ]
    sizes = [(12, 2), (3, 3), (3, 3), (6, 4), (6, 5), (6, 5), (6, 5), (6, 5), (6, 5), (12, 6)]
    layout, x, y, row_h = [], 0, 0, 0
    for chart, (width, height) in zip(charts, sizes):
        if x + width > 12:
            x, y, row_h = 0, y + row_h, 0
        layout.append({"i": chart["id"], "x": x, "y": y, "w": width, "h": height})
        x += width
        row_h = max(row_h, height)
    return {
        "name": f"CRISPR Pipeline — {run_name}",
        "description": f"Live process, resource, tool, and QC status for execution {run_id}.",
        "owner": "X-AXIOM-EVERYONE", "datasets": [dataset], "refreshTime": 60,
        "schemaVersion": 2, "timeWindowStart": "qr-now-24h", "timeWindowEnd": "qr-now",
        "charts": charts, "layout": layout,
    }


def read_new_trace_rows(path: Path, seen: set[tuple[str, str, str]]) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return []
    fresh = []
    for row in rows:
        key = (row.get("task_id", ""), row.get("hash", ""), row.get("attempt", ""))
        if key not in seen and (row.get("status", "").upper() in TERMINAL_STATES):
            seen.add(key)
            fresh.append(row)
    return fresh


def read_new_hook_events(path: Path, offset: int, base: dict[str, Any]) -> tuple[list[dict[str, Any]], int]:
    if not path.exists():
        return [], offset
    try:
        size = path.stat().st_size
        if size < offset:
            offset = 0
        events = []
        with path.open(encoding="utf-8", errors="replace") as handle:
            handle.seek(offset)
            for line in handle:
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue
                events.append({**item, **base, "message": safe_text(item.get("message", ""))})
            offset = handle.tell()
        return events, offset
    except OSError:
        return [], offset


def flatten_qc_metrics(path: Path, base: dict[str, Any]) -> list[dict[str, Any]]:
    if not path.exists() or path.stat().st_size > 10_000_000:
        return []
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    events: list[dict[str, Any]] = []

    def walk(value: Any, key: str = "") -> None:
        if len(events) >= 500:
            return
        if isinstance(value, dict):
            for child_key, child in value.items():
                walk(child, f"{key}.{child_key}".strip("."))
        elif isinstance(value, (int, float, bool)) and not isinstance(value, str):
            event = dict(base)
            event.update({"_time": utc_now(), "event_type": "qc_metric", "status": "available",
                          "metric": key, "value": value, "unit": "", "source": path.name,
                          "message": f"QC metric {key}"})
            events.append(event)

    walk(payload)
    return events


def dependency_events(path: Path, base: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract a bounded process-edge list from Nextflow's Graphviz DAG."""
    if not path.exists() or path.stat().st_size > 5_000_000:
        return []
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    labels = dict(re.findall(r'^\s*([A-Za-z0-9_]+)\s*\[label\s*=\s*"([^"]+)"', content, re.M))
    edges = re.findall(r'^\s*([A-Za-z0-9_]+)\s*->\s*([A-Za-z0-9_]+)', content, re.M)
    events = []
    for upstream, downstream in edges[:1000]:
        source = labels.get(upstream, upstream)
        target = labels.get(downstream, downstream)
        event = dict(base)
        event.update({"_time": utc_now(), "event_type": "process_dependency", "status": "DEFINED",
                      "upstream_process": safe_text(source, 512), "downstream_process": safe_text(target, 512),
                      "message": f"{safe_text(source, 256)} -> {safe_text(target, 256)}"})
        events.append(event)
    return events


def run_command(args: argparse.Namespace) -> int:
    run_id = args.run_id or str(uuid.uuid4())
    run_name = args.run_name or f"nextflow-{run_id[:8]}"
    state_dir = Path(args.state_dir).resolve()
    trace_path = Path(args.trace).resolve()
    dag_path = Path(args.dag).resolve()
    hook_path = Path(args.hook_events).resolve()
    outdir = Path(args.outdir).resolve()
    base = {"service": "crispr_pipeline", "dataset": args.dataset, "run_id": run_id,
            "run_name": run_name, "execution_date": utc_now()[:10],
            "pipeline_revision": args.pipeline_revision or "unknown"}
    client = AxiomClient(args.dataset, args.token_env, args.ingest_url, args.api_url,
                         args.timeout, args.max_bytes, args.max_event_bytes, state_dir)
    env = os.environ.copy()
    env.update({"AXIOM_RUN_ID": run_id, "AXIOM_RUN_NAME": run_name,
                "AXIOM_EVENT_FILE": str(hook_path)})
    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("axiom_telemetry.py run requires a command after --")

    started = time.monotonic()
    process = subprocess.Popen(command, env=env)
    old_handlers: dict[int, Any] = {}
    for signum in (signal.SIGINT, signal.SIGTERM):
        old_handlers[signum] = signal.getsignal(signum)
        signal.signal(signum, lambda sig, _frame, p=process: p.send_signal(sig))

    dashboard_uid = client.create_dashboard(run_id, run_name) if args.create_dashboard else None
    client.ingest([{**base, "_time": utc_now(), "event_type": "run_started", "status": "RUNNING",
                    "message": f"Nextflow run {run_name} started", "dashboard_uid": dashboard_uid or ""}])
    seen: set[tuple[str, str, str]] = set()
    hook_offset = 0
    last_heartbeat = 0.0
    return_code = 1
    try:
        while process.poll() is None:
            rows = read_new_trace_rows(trace_path, seen)
            if rows:
                client.ingest(trace_event(row, base) for row in rows)
            hook_events, hook_offset = read_new_hook_events(hook_path, hook_offset, base)
            if hook_events:
                client.ingest(hook_events)
            now = time.monotonic()
            if now - last_heartbeat >= args.heartbeat_seconds:
                client.ingest([{**base, "_time": utc_now(), "event_type": "heartbeat", "status": "RUNNING",
                                "completed_tasks": len(seen), "elapsed_seconds": int(now - started),
                                "message": f"Run active; {len(seen)} terminal tasks observed"}])
                last_heartbeat = now
            time.sleep(args.poll_seconds)
        return_code = process.wait()
    except Exception as exc:
        client.warn(f"sidecar monitoring failed; waiting for Nextflow unchanged: {exc}")
        return_code = process.wait()
    try:
        rows = read_new_trace_rows(trace_path, seen)
        if rows:
            client.ingest(trace_event(row, base) for row in rows)
        hook_events, hook_offset = read_new_hook_events(hook_path, hook_offset, base)
        if hook_events:
            client.ingest(hook_events)
        qc_path = outdir / "pipeline_qc_metrics.json"
        if not qc_path.exists():
            qc_path = outdir / "pipeline_dashboard" / "pipeline_qc_metrics.json"
        client.ingest(flatten_qc_metrics(qc_path, base))
        client.ingest(dependency_events(dag_path, base))
        client.ingest([{**base, "_time": utc_now(), "event_type": "run_completed",
                        "status": "SUCCEEDED" if return_code == 0 else "FAILED",
                        "exit_code": return_code, "completed_tasks": len(seen),
                        "elapsed_seconds": int(time.monotonic() - started),
                        "telemetry_bytes": client.sent_bytes, "telemetry_dropped_events": client.dropped_events,
                        "message": f"Nextflow run {run_name} exited with code {return_code}"}])
    except Exception as exc:
        client.warn(f"final telemetry flush failed; preserving Nextflow exit code: {exc}")
    finally:
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)
    return return_code


def emit_event(args: argparse.Namespace) -> int:
    base = {"service": "crispr_pipeline", "dataset": args.dataset,
            "run_id": args.run_id or os.environ.get("AXIOM_RUN_ID", "unknown"),
            "run_name": args.run_name or os.environ.get("AXIOM_RUN_NAME", "unknown")}
    client = AxiomClient(args.dataset, args.token_env, args.ingest_url, args.api_url,
                         args.timeout, args.max_bytes, args.max_event_bytes, Path(args.state_dir))
    client.ingest([{**base, "_time": utc_now(), "event_type": args.event_type,
                    "status": args.status, "message": safe_text(args.message)}])
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--dataset", default="crispr_pipeline")
    common.add_argument("--token-env", default="AXIOM_IGVF")
    common.add_argument("--ingest-url", default="https://us-east-1.aws.edge.axiom.co/v1/ingest/{dataset}")
    common.add_argument("--api-url", default="https://api.axiom.co")
    common.add_argument("--timeout", type=float, default=5.0)
    common.add_argument("--max-bytes", type=int, default=DEFAULT_MAX_BYTES)
    common.add_argument("--max-event-bytes", type=int, default=DEFAULT_MAX_EVENT_BYTES)
    common.add_argument("--state-dir", default=".axiom_telemetry")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    run = subparsers.add_parser("run", parents=[common])
    run.add_argument("--run-id", default="")
    run.add_argument("--run-name", default="")
    run.add_argument("--pipeline-revision", default="")
    run.add_argument("--trace", required=True)
    run.add_argument("--dag", required=True)
    run.add_argument("--hook-events", required=True)
    run.add_argument("--outdir", required=True)
    run.add_argument("--poll-seconds", type=float, default=2.0)
    run.add_argument("--heartbeat-seconds", type=float, default=60.0)
    run.add_argument("--create-dashboard", action=argparse.BooleanOptionalAction, default=True)
    run.add_argument("command", nargs=argparse.REMAINDER)
    emit = subparsers.add_parser("emit", parents=[common])
    emit.add_argument("--run-id", default="")
    emit.add_argument("--run-name", default="")
    emit.add_argument("--event-type", required=True)
    emit.add_argument("--status", required=True)
    emit.add_argument("--message", default="")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return run_command(args) if args.subcommand == "run" else emit_event(args)
    except Exception as exc:
        print(f"WARN: Axiom telemetry wrapper failed open: {safe_text(exc)}", file=sys.stderr)
        return 1 if args.subcommand == "run" else 0


if __name__ == "__main__":
    raise SystemExit(main())
