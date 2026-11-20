#!/usr/bin/env python3
"""Carry published pre-inference QC into a separate continuation report."""

import argparse
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

from render_wandb_pipeline_dashboard import CHECKPOINT_FAMILIES, image_family, render

PUBLISHED_QC_DIRS = (
    "measurement_set_qc", "guide_mapping_qc", "guide_assignment_qc",
    "postconcat_embedding_qc", "clone_removal", "hashing_qc", "sequencing_saturation",
)
TABLES = {
    "guide_position_table.csv", "measurement_set_qc_metrics.tsv", "measurement_set_qc_filter_flow.tsv",
    "measurement_filter_flow.tsv", "guide_mapping_qc.json", "gene_metrics.tsv",
    "guide_metrics.tsv", "guide_mapping_qc.tsv", "guide_assignment_filter_flow.tsv",
    "hto_filter_flow.tsv", "clone_metrics.tsv", "sequencing_saturation_metrics.tsv",
    "embedding_qc_metrics.json", "pipeline_qc_metrics.json",
}
MAX_QC_FILE_BYTES = 32_000_000


def copy_metric(source, target):
    if source.is_file() and source.stat().st_size <= MAX_QC_FILE_BYTES:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def copy_qc_dir(source, target):
    if not source.is_dir():
        return
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if path.name in TABLES or (
            path.suffix == ".png"
            and image_family(Path(source.name) / relative) in (*CHECKPOINT_FAMILIES, "final")
        ):
            copy_metric(path, target / relative)


def snapshot(args):
    args.output.mkdir()
    copy_metric(args.manifest, args.output / "checkpoint.json")
    copy_metric(args.metrics, args.output / "pipeline_qc_metrics.json")
    copy_metric(args.seqspec, args.output / "guide_position_table.csv")
    for qc_dir in args.qc_dirs:
        copy_qc_dir(qc_dir, args.output / qc_dir.name)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def restore(args):
    args.output.mkdir()
    report = read_json(args.report)
    checkpoint = args.source / "checkpoint_qc"
    checkpoint_record = read_json(checkpoint / "checkpoint.json")
    if checkpoint_record and (
        checkpoint_record.get("run_name"),
        checkpoint_record.get("mudata_sha256"),
    ) != (report["parent_run_name"], report["parent_sha256"]):
        raise ValueError("the dashboard snapshot does not match the validated checkpoint")

    if checkpoint_record:
        copy_qc_dir(checkpoint, args.output)
    else:
        # Older checkpoints may have pre-inference plots but no compiled metrics.
        for name in PUBLISHED_QC_DIRS:
            copy_qc_dir(args.source / name, args.output / name)
        copy_qc_dir(args.source / "seqspeccheck", args.output / "seqspeccheck")
        seqspec = args.source / "seqspeccheck" / "guide_position_table.csv"
        copy_metric(seqspec, args.output / "guide_position_table.csv")

    copy_metric(args.guide_report, args.output / "guide_metadata.validation.json")
    report["qc_snapshot"] = bool(checkpoint_record)
    report["qc_artifacts"] = any(path.is_file() for path in args.output.rglob("*"))
    (args.output / "continuation.json").write_text(json.dumps(report), encoding="utf-8")


def complete(args):
    import mudata
    from create_dashboard_df import create_inference_blocks

    preview = args.root / "inference_preview.html"
    mdata = mudata.read_h5mu(args.mudata, backed="r")
    try:
        headings = []
        for block in create_inference_blocks(mdata, use_default=True, limit=25):
            label = block["description"].iloc[0]
            table = block["table"].iloc[0].to_html(index=False)
            headings.append(f"<h3>Guide Inference</h3><p>{label}</p>{table}")
        preview.write_text("".join(headings), encoding="utf-8")
    finally:
        mdata.file.close()

    unavailable = args.root / ".not-available"
    document = render(SimpleNamespace(
        trace=unavailable, run_id=args.run_name, run_name=args.run_name, status="completed",
        guide_report=unavailable, seqspec_table=unavailable, seqspec_image=None,
        artifact_dir=args.root, final_dashboard_html=preview,
    ))
    (args.root / "dashboard.html").write_text(document, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    save = subparsers.add_parser("snapshot")
    for name in ("manifest", "metrics", "seqspec", "output"):
        save.add_argument(f"--{name}", required=True, type=Path)
    save.add_argument("qc_dirs", nargs="+", type=Path)
    save.set_defaults(action=snapshot)

    load = subparsers.add_parser("restore")
    for name in ("report", "source", "guide-report", "output"):
        load.add_argument(f"--{name}", required=True, type=Path)
    load.set_defaults(action=restore)

    final = subparsers.add_parser("complete")
    final.add_argument("--root", required=True, type=Path)
    final.add_argument("--mudata", required=True, type=Path)
    final.add_argument("--run-name", required=True)
    final.set_defaults(action=complete)

    args = parser.parse_args()
    args.action(args)


if __name__ == "__main__":
    main()
