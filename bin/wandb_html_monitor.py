#!/usr/bin/env python3
"""Fail-open W&B publisher for one live CRISPR Pipeline HTML dashboard."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

from render_wandb_pipeline_dashboard import render


def warn(message: str) -> None:
    print(f"WARN: W&B HTML telemetry: {message}", file=sys.stderr, flush=True)


def first_file(root: Path, patterns: list[str]) -> Path:
    for pattern in patterns:
        matches = sorted(root.glob(pattern)) if root.exists() else []
        if matches:
            return matches[-1]
    return root / ".not-available"


def discovered_paths(outdir: Path, run_name: str) -> dict[str, Path]:
    return {
        "guide_report": first_file(outdir, [
            f"pipeline_info/{run_name}/guide_metadata.validation.json",
            "pipeline_info/*/guide_metadata.validation.json",
        ]),
        "seqspec_table": first_file(outdir, ["pipeline_outputs/seqspeccheck/guide_position_table.csv"]),
        "seqspec_image": first_file(outdir, ["pipeline_outputs/seqspeccheck/**/*seqSpec*plots.png"]),
        "qc_metrics_json": first_file(outdir, ["pipeline_qc_metrics.json"]),
        # Search the complete published output tree.  This lets the advanced
        # execution dashboard expose QC plots as soon as their producing
        # process publishes them, before the final dashboard is assembled.
        "artifact_dir": outdir,
        "final_dashboard_html": outdir / "pipeline_dashboard" / "dashboard.html",
    }


def input_signature(paths: dict[str, Path], trace: Path, status_file: Path) -> str:
    records = []
    for path in [trace, status_file, *paths.values()]:
        if path.exists():
            stat = path.stat()
            records.append((str(path), stat.st_size, stat.st_mtime_ns))
    return hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()


def read_final_status(status_file: Path) -> tuple[str, bool]:
    if not status_file.exists():
        return "running", False
    try:
        payload = json.loads(status_file.read_text(encoding="utf-8"))
        return str(payload.get("status", "failed")), True
    except (OSError, ValueError) as error:
        warn(f"cannot read status file; using failed: {error}")
        return "failed", True


def render_snapshot(args: argparse.Namespace, status: str, final: bool) -> int:
    paths = discovered_paths(args.outdir, args.run_name)
    render_args = SimpleNamespace(
        trace=args.trace,
        run_id=args.source_run_id,
        run_name=args.run_name,
        status=status,
        guide_report=paths["guide_report"],
        seqspec_table=paths["seqspec_table"],
        seqspec_image=paths["seqspec_image"],
        qc_metrics_json=paths["qc_metrics_json"],
        artifact_dir=paths["artifact_dir"],
        # The final pipeline dashboard is a bounded result-table source for
        # the advanced execution dashboard, never a replacement for it.
        final_dashboard_html=paths["final_dashboard_html"],
        nextflow_log=args.nextflow_log,
        tail_lines=args.tail_lines,
        max_image_bytes=args.max_image_bytes,
    )
    document = render(render_args)
    args.dashboard_html.parent.mkdir(parents=True, exist_ok=True)
    args.dashboard_html.write_text(document, encoding="utf-8")
    return args.dashboard_html.stat().st_size


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="crispr-pipeline")
    parser.add_argument("--entity", default="")
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--wandb-run-id", default="")
    parser.add_argument("--replace-run", default="true")
    parser.add_argument("--publish-live-html", default="true")
    parser.add_argument("--token-env", default="WB_IGVF")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--nextflow-log", type=Path, required=True)
    parser.add_argument("--status-file", type=Path, required=True)
    parser.add_argument("--dashboard-html", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=30)
    parser.add_argument("--max-total-bytes", type=int, default=20_000_000)
    parser.add_argument("--max-final-html-bytes", type=int, default=50_000_000)
    parser.add_argument("--max-image-bytes", type=int, default=10_000_000)
    parser.add_argument("--tail-lines", type=int, default=30)
    args = parser.parse_args()

    token = os.environ.get(args.token_env, "")
    if not token:
        warn(f"{args.token_env} is not set; pipeline continues without W&B")
        return 0
    os.environ["WANDB_API_KEY"] = token
    try:
        import wandb
        publish_live_html = args.publish_live_html.lower() in {"1", "true", "yes"}
        replace_run = args.replace_run.lower() in {"1", "true", "yes"}
        if replace_run and not (args.wandb_run_id and args.entity):
            warn("visible-run replacement requires --wandb-run-id and --entity; using one run")
            replace_run = False
        api = wandb.Api(timeout=30) if replace_run else None
    except Exception as error:
        warn(f"initialization failed; pipeline continues: {error}")
        return 0

    sent_bytes = 0
    step = 0
    previous = ""
    run = None
    old_run_ids: list[str] = []
    if replace_run and api is not None:
        try:
            for candidate in api.runs(f"{args.entity}/{args.project}"):
                series = candidate.config.get("dashboard_series_id", "")
                if candidate.id == args.wandb_run_id or series == args.wandb_run_id:
                    old_run_ids.append(candidate.id)
        except Exception as error:
            warn(f"cannot inventory prior dashboard runs: {error}")

    def open_run(run_id: str, resume: str):
        return wandb.init(
            project=args.project,
            entity=args.entity or None,
            id=run_id or None,
            resume=resume if run_id else None,
            name=args.run_name,
            job_type="pipeline-execution-dashboard",
            tags=["crispr-pipeline", "html-dashboard", "live"],
            config={
                "telemetry_layout": "single-visible-html",
                "source_run_id": args.source_run_id,
                "dashboard_series_id": args.wandb_run_id,
            },
            settings=wandb.Settings(init_timeout=20),
        )

    try:
        if not replace_run:
            run = open_run(args.wandb_run_id, "allow")
        while True:
            status, final = read_final_status(args.status_file)
            paths = discovered_paths(args.outdir, args.run_name)
            signature = input_signature(paths, args.trace, args.status_file)
            if signature != previous or final:
                try:
                    if final or publish_live_html:
                        size = render_snapshot(args, status, final)
                        within_budget = size <= args.max_final_html_bytes
                        if within_budget:
                            if replace_run:
                                # W&B only materializes a visible HTML panel
                                # from history. Each refresh therefore gets a
                                # fresh one-point run; preceding members of the
                                # dashboard series are removed after success.
                                suffix = f"{time.time_ns():x}"[-14:]
                                visible_id = f"{args.wandb_run_id[:70]}-{suffix}"
                                run = open_run(visible_id, "never")
                            assert run is not None
                            run.log({
                                "pipeline/main_execution": wandb.Html(
                                    str(args.dashboard_html), inject=False
                                )
                            }, step=0)
                            run.summary.update({
                                "pipeline_status": status,
                                "dashboard_status": (
                                    "FULL_QC_PUBLISHED" if final else "LIVE"
                                ),
                                "dashboard_updates": step + 1,
                                "dashboard_uploaded_bytes": sent_bytes + size,
                                "source_run_id": args.source_run_id,
                            })
                            if replace_run:
                                run.finish()
                                run = None
                                for old_id in old_run_ids:
                                    if old_id == visible_id:
                                        continue
                                    try:
                                        api.run(
                                            f"{args.entity}/{args.project}/{old_id}"
                                        ).delete()
                                    except Exception as error:
                                        warn(f"cannot remove prior dashboard {old_id}: {error}")
                                old_run_ids = [visible_id]
                            sent_bytes += size
                            step += 1
                        else:
                            warn(
                                f"upload limit reached ({sent_bytes} bytes sent); "
                                f"skipping {size}-byte snapshot"
                            )
                    previous = signature
                except Exception as error:
                    warn(f"render/upload failed; pipeline continues: {error}")
            if final:
                break
            time.sleep(max(args.poll_seconds, 1))
        if run is not None:
            run.summary.update({
                "dashboard_status": "FULL_QC_PUBLISHED" if status == "completed" else "FAILED",
                "dashboard_updates": step,
                "dashboard_uploaded_bytes": sent_bytes,
                "source_run_id": args.source_run_id,
            })
    except Exception as error:
        warn(f"monitor failed; pipeline continues: {error}")
    finally:
        try:
            if run is not None:
                run.finish()
        except Exception as error:
            warn(f"finish failed; pipeline continues: {error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
