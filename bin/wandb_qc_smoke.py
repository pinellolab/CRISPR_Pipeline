#!/usr/bin/env python3
"""Publish one self-contained CRISPR Pipeline execution dashboard to W&B."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="crispr-pipeline")
    parser.add_argument("--entity", default="")
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--token-env", default="WB_IGVF")
    parser.add_argument("--dashboard-html", type=Path, required=True)
    parser.add_argument("--source-run-id", default="")
    parser.add_argument("--max-bytes", type=int, default=20_000_000)
    args = parser.parse_args()

    token = os.environ.get(args.token_env, "")
    if not token:
        raise SystemExit(f"{args.token_env} is not set")
    os.environ["WANDB_API_KEY"] = token

    import wandb

    if not args.dashboard_html.is_file():
        raise SystemExit(f"dashboard HTML does not exist: {args.dashboard_html}")
    dashboard_bytes = args.dashboard_html.stat().st_size
    if dashboard_bytes > args.max_bytes:
        raise SystemExit(f"dashboard HTML exceeds --max-bytes ({args.max_bytes})")

    run = wandb.init(
        project=args.project,
        entity=args.entity or None,
        name=args.run_name,
        job_type="pipeline-execution-dashboard",
        tags=["crispr-pipeline", "html-dashboard"],
        config={"telemetry_layout": "single-html", "source_run_id": args.source_run_id},
    )
    run.log({
        "pipeline/main_execution": wandb.Html(str(args.dashboard_html), inject=False),
    }, step=0)
    run.summary["dashboard_status"] = "PUBLISHED"
    run.summary["dashboard_bytes"] = dashboard_bytes
    run.summary["source_run_id"] = args.source_run_id
    print(json.dumps({"run_id": run.id, "run_url": run.url, "entity": run.entity,
                      "project": run.project, "dashboard_bytes": dashboard_bytes}))
    run.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
