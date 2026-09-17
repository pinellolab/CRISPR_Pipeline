#!/usr/bin/env python3
"""Fail-open W&B publisher for one live CRISPR Pipeline HTML dashboard."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import mimetypes
import os
import re
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
        "artifact_dir": outdir / "pipeline_dashboard",
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


def inline_dashboard_assets(dashboard_html: Path) -> str:
    """Return the canonical white dashboard with local image assets embedded."""
    source = dashboard_html.read_text(encoding="utf-8", errors="replace")
    root = dashboard_html.parent.resolve()
    cache: dict[Path, str] = {}

    def replace(match: re.Match[str]) -> str:
        attribute, quote, value = match.groups()
        if value.startswith(("data:", "http:", "https:", "#", "javascript:")):
            return match.group(0)
        candidate = (root / value).resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            return match.group(0)
        if not candidate.is_file():
            return match.group(0)
        if candidate not in cache:
            media_type = mimetypes.guess_type(candidate.name)[0] or "application/octet-stream"
            encoded = base64.b64encode(candidate.read_bytes()).decode("ascii")
            cache[candidate] = f"data:{media_type};base64,{encoded}"
        return f"{attribute}={quote}{cache[candidate]}{quote}"

    # Figures are loaded lazily from data-imgsrc; navigation icons use src.
    # Embed both so W&B receives the complete dashboard as one HTML object.
    return re.sub(r"\b(src|data-imgsrc)=([\"'])([^\"']+)\2", replace, source)


def render_snapshot(args: argparse.Namespace, status: str, final: bool) -> int:
    paths = discovered_paths(args.outdir, args.run_name)
    if final and paths["final_dashboard_html"].is_file():
        document = inline_dashboard_assets(paths["final_dashboard_html"])
        args.dashboard_html.parent.mkdir(parents=True, exist_ok=True)
        args.dashboard_html.write_text(document, encoding="utf-8")
        return args.dashboard_html.stat().st_size
    render_args = SimpleNamespace(
        trace=args.trace,
        run_id=args.source_run_id,
        run_name=args.run_name,
        status=status,
        guide_report=paths["guide_report"],
        seqspec_table=paths["seqspec_table"],
        seqspec_image=paths["seqspec_image"] if final else None,
        qc_metrics_json=paths["qc_metrics_json"],
        artifact_dir=paths["artifact_dir"] if final else None,
        final_dashboard_html=paths["final_dashboard_html"] if final else None,
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
    parser.add_argument("--replace-run", default="false")
    parser.add_argument("--publish-live-html", default="false")
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
        replace_run = args.replace_run.lower() in {"1", "true", "yes"}
        publish_live_html = args.publish_live_html.lower() in {"1", "true", "yes"}
        if replace_run and args.wandb_run_id and args.entity:
            try:
                wandb.Api(timeout=30).run(
                    f"{args.entity}/{args.project}/{args.wandb_run_id}"
                ).delete()
            except wandb.errors.CommError:
                pass
        init_options = {}
        if args.wandb_run_id:
            init_options.update(id=args.wandb_run_id, resume="allow")
        run = wandb.init(
            project=args.project,
            entity=args.entity or None,
            name=args.run_name,
            job_type="pipeline-execution-dashboard",
            tags=["crispr-pipeline", "html-dashboard", "live"],
            config={"telemetry_layout": "single-html", "source_run_id": args.source_run_id},
            settings=wandb.Settings(init_timeout=20),
            **init_options,
        )
    except Exception as error:
        warn(f"initialization failed; pipeline continues: {error}")
        return 0

    sent_bytes = 0
    step = 0
    previous = ""
    try:
        while True:
            status, final = read_final_status(args.status_file)
            paths = discovered_paths(args.outdir, args.run_name)
            signature = input_signature(paths, args.trace, args.status_file)
            if signature != previous or final:
                try:
                    run.summary["pipeline_status"] = status
                    run.summary["dashboard_status"] = (
                        "BUILDING_FULL_QC" if not final else "FINALIZING"
                    )
                    if final or publish_live_html:
                        size = render_snapshot(args, status, final)
                        full_dashboard = final and paths["final_dashboard_html"].is_file()
                        within_budget = (
                            size <= args.max_final_html_bytes
                            if full_dashboard
                            else sent_bytes + size <= args.max_total_bytes
                        )
                        if within_budget:
                            run.log({"pipeline/main_execution": wandb.Html(str(args.dashboard_html), inject=False)}, step=step)
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
        run.summary["dashboard_status"] = "FULL_QC_PUBLISHED" if status == "completed" else "FAILED"
        run.summary["dashboard_updates"] = step
        run.summary["dashboard_uploaded_bytes"] = sent_bytes
        run.summary["source_run_id"] = args.source_run_id
    except Exception as error:
        warn(f"monitor failed; pipeline continues: {error}")
    finally:
        try:
            run.finish()
        except Exception as error:
            warn(f"finish failed; pipeline continues: {error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
