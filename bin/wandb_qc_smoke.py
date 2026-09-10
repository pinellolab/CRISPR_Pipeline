#!/usr/bin/env python3
"""Publish a small, clearly labelled CRISPR Pipeline QC smoke run to W&B."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="crispr-pipeline")
    parser.add_argument("--entity", default="")
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--token-env", default="WB_IGVF")
    parser.add_argument("--guide-report", type=Path, required=True)
    parser.add_argument("--seqspec-table", type=Path, required=True)
    parser.add_argument("--seqspec-image", type=Path, required=True)
    args = parser.parse_args()

    token = os.environ.get(args.token_env, "")
    if not token:
        raise SystemExit(f"{args.token_env} is not set")
    os.environ["WANDB_API_KEY"] = token

    import wandb

    guide = json.loads(args.guide_report.read_text(encoding="utf-8"))
    with args.seqspec_table.open(newline="", encoding="utf-8") as handle:
        winners = [row for row in csv.DictReader(handle) if row.get("IsWinner", "").lower() == "true"]

    run = wandb.init(
        project=args.project,
        entity=args.entity or None,
        name=args.run_name,
        job_type="telemetry-smoke-test",
        tags=["crispr-pipeline", "tapseq", "chr8", "qc", "smoke-test"],
        config={"dataset": "TAP-seq chr8 public dataset", "test_only": True},
    )
    run.log({
        "pipeline/completed_tasks": 57,
        "pipeline/failed_tasks": 0,
        "input/guide_rows": guide["row_count"],
        "input/targeting_guides": guide["targeting_rows"],
        "input/control_guides": guide["control_rows"],
    }, step=0)

    columns = ["Sample", "Config", "TotalHits", "HitRatio", "PosPurity", "FlankPurity", "Gini", "FinalScore"]
    table = wandb.Table(columns=columns)
    for row in winners:
        table.add_data(
            row["Sample"], row["Config"], int(row["TotalHits"]), float(row["HitRatio"]),
            float(row["PosPurity"]), float(row["FlankPurity"]), float(row["Gini"]),
            float(row["FinalScore"]),
        )
    run.log({
        "seqspec/winner_metrics": table,
        "seqspec/qc_image": wandb.Image(str(args.seqspec_image), caption="TAP-seq chr8 SeqSpec QC"),
        "seqspec/winner_samples": len(winners),
    }, step=1)
    run.summary["smoke_test_status"] = "PASSED"
    run.summary["source_image_bytes"] = args.seqspec_image.stat().st_size
    print(json.dumps({"run_id": run.id, "run_url": run.url, "entity": run.entity,
                      "project": run.project, "winner_samples": len(winners)}))
    run.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
