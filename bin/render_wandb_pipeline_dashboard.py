#!/usr/bin/env python3
"""Render a self-contained, clickable CRISPR Pipeline execution dashboard."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


FAMILIES = [
    ("input", "Input QC", "Samples, guides and provenance"),
    ("seqspec", "SeqSpec", "Read structure and capture QC"),
    ("mapping", "Mapping", "RNA and guide quantification"),
    ("preprocessing", "Preprocessing", "Cell and gene filtering"),
    ("mudata", "MuData", "Modalities assembled"),
    ("guide_assignment", "Guide assignment", "Guide-to-cell calls"),
    ("inference", "Inference", "SCEPTRE and Perturbo"),
    ("evaluation", "Evaluation", "Controls and benchmarking"),
    ("final", "Final dashboard", "Published report and artifacts"),
]


def family_for(process: str) -> str:
    value = process.lower()
    leaf = value.split(":")[-1]
    if "seqspec" in value:
        return "seqspec"
    if any(key in value for key in ("mapping_rna_pipeline", "mapping_guide_pipeline", "mapping_hashing_pipeline")):
        if any(key in leaf for key in ("downloadreference", "seqspecparser", "createguideref", "createhashingref")):
            return "input"
        return "mapping"
    if any(key in value for key in ("preprocessing_pipeline", "preprocessanndata", "doublets", "filter_hashing")):
        return "preprocessing"
    if any(key in value for key in ("createmudata", "anndata_concat", "mudata_concat", "hashing_concat")):
        return "mudata"
    if "guide_assignment" in value or "prepare_assignment" in value:
        return "guide_assignment"
    if any(key in value for key in ("inference", "sceptre_chunk", "perturbo", "mergedresults", "catalog", "mergemudata")):
        return "inference"
    if any(key in value for key in ("evaluation", "additional_qc", "benchmark")):
        return "evaluation"
    if "dashboard" in value or "publishfiles" in value:
        return "final"
    return "input"


def read_trace(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: Any) -> float | None:
    try:
        return float(str(value).replace("%", "").replace(",", ""))
    except (TypeError, ValueError):
        return None


def duration_seconds(value: str) -> float:
    import re
    factors = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}
    return sum(float(n) * factors[u] for n, u in re.findall(r"([0-9.]+)\s*(ms|s|m|h|d)", value or ""))


def family_state(rows: list[dict[str, str]], family: str, run_status: str) -> dict[str, Any]:
    selected = [row for row in rows if family_for(row.get("process", "")) == family]
    statuses = Counter(row.get("status", "UNKNOWN").upper() for row in selected)
    failed = statuses["FAILED"] + statuses["ABORTED"]
    if failed:
        status = "failed"
    elif selected:
        status = "completed"
    else:
        status = "pending"
    if run_status.lower() == "running" and selected and family == next(
        (fid for fid, _, _ in reversed(FAMILIES) if any(family_for(r.get("process", "")) == fid for r in rows)), "input"
    ):
        status = "running"
    return {
        "status": status,
        "rows": selected,
        "completed": statuses["COMPLETED"] + statuses["CACHED"],
        "failed": failed,
        "cached": statuses["CACHED"],
        "runtime": sum(duration_seconds(row.get("realtime") or row.get("duration", "")) for row in selected),
    }


def fmt_seconds(seconds: float) -> str:
    if seconds >= 3600:
        return f"{seconds / 3600:.1f} h"
    if seconds >= 60:
        return f"{seconds / 60:.1f} min"
    return f"{seconds:.1f} s"


def metric_card(label: str, value: Any, detail: str = "") -> str:
    return (
        '<div class="metric"><div class="metric-label">' + html.escape(label) + '</div>'
        '<div class="metric-value">' + html.escape(str(value)) + '</div>'
        '<div class="metric-detail">' + html.escape(detail) + '</div></div>'
    )


def display_value(value: Any) -> str:
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if abs(value) >= 1000:
            return f"{value:,.0f}"
        if abs(value) >= 10:
            return f"{value:,.1f}"
        return f"{value:.3f}"
    return str(value)


def data_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]], limit: int = 100) -> str:
    if not rows:
        return '<div class="empty">Metrics are not available yet.</div>'
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = []
    for row in rows[:limit]:
        body.append("<tr>" + "".join(
            f"<td>{html.escape(display_value(row.get(key)))}</td>" for key, _ in columns
        ) + "</tr>")
    return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def overall_row(section: dict[str, Any]) -> dict[str, Any]:
    rows = section.get("rows", []) if isinstance(section, dict) else []
    return next((row for row in rows if row.get("batch") == "all"), rows[0] if rows else {})


def qc_metrics_content(data: dict[str, Any], family: str) -> str:
    if not data:
        return ""
    observed = data.get("observed_metrics", {})
    additional = observed.get("additional_qc", {})
    summary = data.get("dashboard_summary", {}).get("filtering_summary", {})

    if family == "input":
        params = data.get("parameters", {}).get("selected_qc_and_inference_params", {})
        rows = [{"parameter": key, "value": value} for key, value in sorted(params.items())]
        return "<h3>Resolved analysis settings</h3>" + data_table(rows, [("parameter", "Parameter"), ("value", "Value")])

    if family == "mapping":
        rows = []
        for item in observed.get("mapping_json", []):
            metrics = item.get("metrics", {})
            rows.append({
                "measurement_set": item.get("measurement_set"), "modality": item.get("modality"),
                "processed": metrics.get("n_processed"), "pseudoaligned": metrics.get("p_pseudoaligned"),
                "unique": metrics.get("p_unique"), "reads_onlist": metrics.get("percentageReadsOnOnlist"),
                "barcodes_onlist": metrics.get("percentageBarcodesOnOnlist"),
                "median_umis": metrics.get("medianUMIsPerBarcode"),
            })
        return "<h3>Mapping and barcode capture</h3>" + data_table(rows, [
            ("measurement_set", "Measurement set"), ("modality", "Modality"),
            ("processed", "Reads processed"), ("pseudoaligned", "Pseudoaligned %"),
            ("unique", "Unique %"), ("reads_onlist", "Reads on-list %"),
            ("barcodes_onlist", "Barcodes on-list %"), ("median_umis", "Median UMI/barcode"),
        ])

    if family == "preprocessing":
        barcode = summary.get("barcode_filter", {})
        diagnostics = summary.get("diagnostics", {})
        gene = additional.get("gene", {})
        overall = overall_row(gene)
        cards = [
            metric_card("Raw RNA barcodes", display_value(summary.get("raw_scRNA_barcodes_total")), "before concatenation"),
            metric_card("Concatenated cells", display_value(summary.get("concatenated_scRNA_cells")), "RNA modality"),
            metric_card("Cells after filter", display_value(barcode.get("cells_after_filter")), barcode.get("method", "barcode filter")),
            metric_card("Median RNA UMI", display_value(overall.get("umi_median")), "final cells"),
            metric_card("Mean RNA UMI", display_value(overall.get("umi_mean")), "final cells"),
            metric_card("Median genes", display_value(diagnostics.get("median_detected_genes_final_cells")), "final cells"),
            metric_card("Mean mito %", display_value(overall.get("mito_mean")), "final cells"),
            metric_card("Mito cutoff %", display_value(summary.get("mitochondrial_filter", {}).get("pct_mito_max")), "configured"),
        ]
        return '<h3>Cell and RNA quality</h3><div class="metrics">' + "".join(cards) + "</div>" + data_table(
            gene.get("rows", []), [("batch", "Batch"), ("n_cells", "Cells"),
            ("umi_median", "Median UMI"), ("umi_mean", "Mean UMI"),
            ("mito_median", "Median mito %"), ("mito_mean", "Mean mito %")])

    if family == "mudata":
        final = summary.get("final_mudata", {})
        intersection = summary.get("modality_intersection", {})
        cards = [
            metric_card("Final cells", display_value(final.get("cells")), "MuData"),
            metric_card("Gene features", display_value(final.get("gene_features")), "RNA modality"),
            metric_card("Guide features", display_value(final.get("guide_features")), "guide modality"),
            metric_card("RNA-guide overlap", display_value(intersection.get("rna_guide_intersection_cells")), "cells"),
            metric_card("Unfiltered overlap", display_value(intersection.get("guide_unfiltered_intersection_fraction")), "fraction"),
            metric_card("Guide assignments", display_value(final.get("total_sgrna_assignment_values")), "total values"),
        ]
        return '<h3>Final multimodal object</h3><div class="metrics">' + "".join(cards) + "</div>"

    if family == "guide_assignment":
        guide = additional.get("guide", {})
        overall = overall_row(guide)
        cards = [
            metric_card("Cells with guide", display_value(overall.get("n_cells_with_guide")), "assigned"),
            metric_card("Assignment rate", display_value(overall.get("frac_cells_with_guide")), "fraction"),
            metric_card("Exactly one guide", display_value(overall.get("n_cells_exactly_1_guide")), "cells"),
            metric_card("Median guide UMI", display_value(overall.get("guide_umi_median")), "per cell"),
            metric_card("Guides per cell", display_value(overall.get("guides_per_cell_mean")), "mean"),
            metric_card("Cells per guide", display_value(overall.get("cells_per_guide_median")), "median"),
            metric_card("Recovered guides", display_value(overall.get("n_guides_total")), "features"),
        ]
        return '<h3>Guide assignment quality</h3><div class="metrics">' + "".join(cards) + "</div>" + data_table(
            guide.get("rows", []), [("batch", "Batch"), ("n_cells", "Cells"),
            ("guide_umi_median", "Median guide UMI"), ("frac_cells_with_guide", "Assigned fraction"),
            ("n_cells_exactly_1_guide", "Exactly one guide"), ("guides_per_cell_mean", "Mean guides/cell")])

    if family == "evaluation":
        intended = overall_row(additional.get("intended_target", {}))
        global_qc = overall_row(additional.get("global_analysis", {}))
        cards = [
            metric_card("Intended tested", display_value(intended.get("n_guides_tested")), "guides"),
            metric_card("Strong knockdowns", display_value(intended.get("n_strong_knockdowns")), "guides"),
            metric_card("Intended significant", display_value(intended.get("n_significant")), "guides"),
            metric_card("Median intended log2FC", display_value(intended.get("median_log2fc")), "effect"),
            metric_card("Global guides tested", display_value(global_qc.get("n_guides_tested")), "guides"),
            metric_card("Significant trans tests", display_value(global_qc.get("total_significant_tests")), "tests"),
            metric_card("Global AUROC", display_value(global_qc.get("auroc")), "validated links"),
            metric_card("Global AUPRC", display_value(global_qc.get("auprc")), "validated links"),
        ]
        rows = [{"analysis": "Intended target", **intended}, {"analysis": "Global", **global_qc}]
        return '<h3>Effect and control evaluation</h3><div class="metrics">' + "".join(cards) + "</div>" + data_table(
            rows, [("analysis", "Analysis"), ("n_guides_tested", "Guides tested"),
            ("n_significant", "Significant"), ("median_log2fc", "Median log2FC"),
            ("auroc", "AUROC"), ("auprc", "AUPRC")])

    if family == "final":
        sources = data.get("sources", {})
        rows = [{"source": key, "value": value} for key, value in sorted(sources.items())]
        return f'<h3>Metric manifest</h3><p>Schema {html.escape(str(data.get("schema_version", "—")))}</p>' + data_table(rows, [("source", "Source"), ("value", "Value")])
    return ""


def image_family(path: Path) -> str:
    value = str(path).lower()
    if "seqspec" in value:
        return "seqspec"
    if any(term in value for term in ("guide_", "guides_", "sgrna", "cells_per_guide", "guides_per_cell")):
        return "guide_assignment"
    if any(term in value for term in ("intended_target", "global_analysis", "evaluation", "volcano")):
        return "evaluation"
    if any(term in value for term in ("scrna", "gene_", "knee_plot")):
        return "preprocessing"
    if any(term in value for term in ("loss_curve", "perturbo", "sceptre")):
        return "inference"
    return "final"


def collect_images(root: Path | None, max_bytes: int = 10_000_000) -> dict[str, list[tuple[Path, str]]]:
    selected: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    if not root or not root.exists():
        return selected
    seen: set[str] = set()
    used = 0
    for path in sorted(root.rglob("*.png")):
        size = path.stat().st_size
        if size > 3_000_000 or used + size > max_bytes:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        used += size
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        selected[image_family(path)].append((path, encoded))
    return selected


def image_gallery(images: list[tuple[Path, str]]) -> str:
    if not images:
        return ""
    figures = "".join(
        f'<figure><img loading="lazy" src="data:image/png;base64,{encoded}" alt="{html.escape(path.stem)}">'
        f'<figcaption>{html.escape(path.stem.replace("_", " "))}</figcaption></figure>'
        for path, encoded in images
    )
    return f'<h3>QC visualizations</h3><div class="gallery">{figures}</div>'


def sanitized_tail(path: Path, lines: int) -> str:
    if not path.exists() or not path.is_file():
        return ""
    raw = path.read_bytes()[-128_000:].decode("utf-8", "replace")
    text = "\n".join(raw.splitlines()[-lines:])
    text = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", text)
    text = re.sub(r"(?i)(authorization\s*[:=]\s*bearer\s+)[^\s]+", r"\1<redacted>", text)
    text = re.sub(r"(?i)\b((?:api_?)?(?:key|token|secret|password))\s*[:=]\s*[^\s]+", r"\1=<redacted>", text)
    text = re.sub(r"\b(?:xaat-|wandb-|axiom-)[A-Za-z0-9._-]{12,}\b", "<redacted-token>", text)
    return text


def failure_content(rows: list[dict[str, str]], nextflow_log: Path | None, tail_lines: int) -> str:
    failed = [row for row in rows if row.get("status", "").upper() in {"FAILED", "ABORTED"}]
    log_tail = sanitized_tail(nextflow_log, max(tail_lines * 2, 40)) if nextflow_log else ""
    if not failed and not any(marker in log_tail for marker in ("ERROR", "Exception", "Error executing")):
        return ""
    blocks = []
    for row in failed:
        workdir = Path(row.get("workdir", ""))
        exitcode = sanitized_tail(workdir / ".exitcode", 1) if workdir.is_dir() else ""
        evidence = []
        for filename in (".command.err", ".command.out"):
            tail = sanitized_tail(workdir / filename, tail_lines) if workdir.is_dir() else ""
            if tail:
                evidence.append(f'<h4>{html.escape(filename)} tail</h4><pre>{html.escape(tail)}</pre>')
        note = ""
        if row.get("status", "").upper() in {"FAILED", "ABORTED"} and exitcode == "0":
            note = '<p class="evidence-note">Task wrapper exit code is 0; the trace failure is consistent with external interruption rather than a tool error.</p>'
        blocks.append(
            '<details open><summary>' + html.escape(row.get("name") or row.get("process", "unknown process")) + '</summary>'
            f'<p>Trace status: <strong>{html.escape(row.get("status", "UNKNOWN"))}</strong> · trace exit: {html.escape(row.get("exit", "—"))} · wrapper exit: {html.escape(exitcode or "—")}</p>'
            + note + "".join(evidence) + "</details>"
        )
    if log_tail:
        blocks.insert(0, f'<details open><summary>Nextflow error / shutdown tail</summary><pre>{html.escape(log_tail)}</pre></details>')
    return '<section class="failure-evidence"><span class="eyebrow">Bounded diagnostic evidence</span><h2>Failure evidence</h2>' + "".join(blocks) + "</section>"


def process_table(rows: list[dict[str, str]]) -> str:
    if not rows:
        return '<div class="empty">No processes observed for this family yet.</div>'
    body = []
    for row in rows[-80:]:
        status = row.get("status", "UNKNOWN").lower()
        body.append(
            "<tr><td>" + html.escape(row.get("process", "").split(":")[-1]) + "</td>"
            "<td>" + html.escape(row.get("name", "")) + "</td>"
            f'<td><span class="pill {status}">{html.escape(row.get("status", "UNKNOWN"))}</span></td>'
            "<td>" + html.escape(row.get("realtime") or row.get("duration", "—")) + "</td>"
            "<td>" + html.escape(row.get("peak_rss", "—")) + "</td></tr>"
        )
    return (
        '<div class="table-wrap"><table><thead><tr><th>Process</th><th>Task</th><th>Status</th>'
        '<th>Runtime</th><th>Peak memory</th></tr></thead><tbody>' + "".join(body) + "</tbody></table></div>"
    )


def seqspec_content(table_path: Path, image_path: Path | None) -> str:
    winners: list[dict[str, str]] = []
    if table_path.exists():
        with table_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            winners = [row for row in csv.DictReader(handle) if row.get("IsWinner", "").lower() == "true"]
    table = '<div class="empty">SeqSpec metrics are not available yet.</div>'
    if winners:
        rows = "".join(
            f"<tr><td>{html.escape(row['Sample'])}</td><td>{html.escape(row['Config'])}</td>"
            f"<td>{html.escape(row['TotalHits'])}</td><td>{float(row['HitRatio']):.3f}</td>"
            f"<td>{float(row['PosPurity']):.3f}</td><td>{float(row['FlankPurity']):.3f}</td>"
            f"<td>{float(row['FinalScore']):.2f}</td></tr>" for row in winners
        )
        table = ('<div class="table-wrap"><table><thead><tr><th>Sample</th><th>Configuration</th>'
                 '<th>Total hits</th><th>Hit ratio</th><th>Position purity</th><th>Flank purity</th>'
                 f'<th>Score</th></tr></thead><tbody>{rows}</tbody></table></div>')
    image = ""
    if image_path and image_path.exists() and image_path.stat().st_size <= 4_000_000:
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        image = f'<figure><img src="data:image/png;base64,{encoded}" alt="SeqSpec QC"><figcaption>SeqSpec read-structure QC</figcaption></figure>'
    return table + image


def render(args: argparse.Namespace) -> str:
    trace_rows = read_trace(args.trace)
    family_states = {family: family_state(trace_rows, family, args.status) for family, _, _ in FAMILIES}
    counts = Counter(row.get("status", "UNKNOWN").upper() for row in trace_rows)
    total_runtime = sum(duration_seconds(row.get("realtime") or row.get("duration", "")) for row in trace_rows)
    guide: dict[str, Any] = {}
    if args.guide_report.exists():
        guide = json.loads(args.guide_report.read_text(encoding="utf-8"))
    qc_data: dict[str, Any] = {}
    qc_metrics_json = getattr(args, "qc_metrics_json", None)
    if qc_metrics_json and qc_metrics_json.exists():
        qc_data = json.loads(qc_metrics_json.read_text(encoding="utf-8"))
    artifact_images = collect_images(
        getattr(args, "artifact_dir", None), getattr(args, "max_image_bytes", 10_000_000)
    )

    graph_nodes = []
    for index, (family, title, subtitle) in enumerate(FAMILIES, start=1):
        state = family_states[family]
        graph_nodes.append(
            f'<button class="node {state["status"]}" data-family="{family}" onclick="selectFamily(\'{family}\')">'
            f'<span class="node-index">{index:02d}</span><span class="node-status"></span>'
            f'<strong>{html.escape(title)}</strong><small>{html.escape(subtitle)}</small>'
            f'<span class="node-count">{state["completed"]} complete · {state["failed"]} failed</span></button>'
        )

    sections = []
    for family, title, subtitle in FAMILIES:
        state = family_states[family]
        cards = [
            metric_card("Completed", state["completed"], "processes"),
            metric_card("Failed", state["failed"], "processes"),
            metric_card("Cached", state["cached"], "processes"),
            metric_card("Task runtime", fmt_seconds(state["runtime"]), "aggregate"),
        ]
        extra = qc_metrics_content(qc_data, family)
        if family == "input" and guide:
            cards.extend([
                metric_card("Guides", guide.get("row_count", "—"), "validated"),
                metric_card("Targeting", guide.get("targeting_rows", "—"), "guides"),
                metric_card("Controls", guide.get("control_rows", "—"), "guides"),
            ])
        if family == "seqspec":
            seqspec_image = None if getattr(args, "artifact_dir", None) else args.seqspec_image
            extra += seqspec_content(args.seqspec_table, seqspec_image)
        extra += image_gallery(artifact_images.get(family, []))
        sections.append(
            f'<section id="family-{family}" class="family-panel"><div class="family-heading">'
            f'<div><span class="eyebrow">Pipeline family</span><h2>{html.escape(title)}</h2>'
            f'<p>{html.escape(subtitle)}</p></div><span class="status-badge {state["status"]}">{state["status"]}</span></div>'
            f'<div class="metrics">{"".join(cards)}</div>{extra}<h3>Processes</h3>{process_table(state["rows"])}</section>'
        )

    current = next((family for family, _, _ in reversed(FAMILIES) if family_states[family]["status"] in {"failed", "running", "completed"}), "input")
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    failures = failure_content(
        trace_rows, getattr(args, "nextflow_log", None), getattr(args, "tail_lines", 30)
    )
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CRISPR Pipeline · {html.escape(args.run_name)}</title>
<style>
:root{{--bg:#07111f;--panel:#0d1b2d;--panel2:#11243a;--line:#29405d;--text:#edf5ff;--muted:#91a7c0;--cyan:#46d9ff;--green:#36d399;--red:#ff647c;--amber:#ffbd59;--grey:#63758a}}
*{{box-sizing:border-box}} body{{margin:0;background:radial-gradient(circle at 20% 0,#12304b 0,transparent 38%),var(--bg);color:var(--text);font:14px/1.5 Inter,ui-sans-serif,system-ui,sans-serif}}
.shell{{max-width:1500px;margin:auto;padding:28px}} header{{display:flex;justify-content:space-between;gap:24px;align-items:flex-start;margin-bottom:20px}} h1{{font-size:27px;margin:4px 0}} h2{{margin:3px 0 0;font-size:24px}} h3{{margin-top:28px}} p{{color:var(--muted);margin:4px 0}} .eyebrow{{color:var(--cyan);font:600 11px ui-monospace,monospace;letter-spacing:.14em;text-transform:uppercase}}
.run-id{{font-family:ui-monospace,monospace;color:var(--muted)}} .live{{display:flex;align-items:center;gap:8px;background:#102840;border:1px solid #28587a;border-radius:99px;padding:8px 12px}} .live i{{width:9px;height:9px;border-radius:50%;background:var(--amber);box-shadow:0 0 12px var(--amber)}} .live.completed i{{background:var(--green);box-shadow:0 0 12px var(--green)}} .live.failed i,.live.interrupted i{{background:var(--red);box-shadow:0 0 12px var(--red)}} .live.running i{{background:var(--cyan);box-shadow:0 0 12px var(--cyan)}}
.summary,.metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:10px;margin:18px 0}} .metric{{background:linear-gradient(145deg,var(--panel2),var(--panel));border:1px solid var(--line);border-radius:12px;padding:14px}} .metric-label{{color:var(--muted);font-size:12px}} .metric-value{{font-size:24px;font-weight:750;margin-top:4px}} .metric-detail{{color:#66809d;font-size:11px}}
.graph-card,.family-panel{{background:rgba(13,27,45,.92);border:1px solid var(--line);border-radius:16px;padding:18px;margin-top:14px;box-shadow:0 18px 55px #0004}} .graph-head{{display:flex;justify-content:space-between;align-items:center}} .graph{{display:flex;align-items:stretch;overflow-x:auto;padding:20px 2px 10px}} .node{{position:relative;flex:0 0 145px;min-height:132px;text-align:left;color:var(--text);background:#102036;border:1px solid var(--line);border-radius:12px;padding:14px;cursor:pointer;transition:.18s}} .node:hover,.node.active{{transform:translateY(-3px);border-color:var(--cyan);box-shadow:0 0 0 2px #46d9ff22}} .node:not(:last-child){{margin-right:29px}} .node:not(:last-child):after{{content:'→';position:absolute;right:-23px;top:49px;color:#55708d;font-size:22px}} .node strong,.node small,.node-count{{display:block}} .node strong{{margin-top:18px}} .node small{{color:var(--muted);font-size:11px;min-height:34px}} .node-count{{font-size:10px;color:#7890aa;margin-top:7px}} .node-index{{font:600 10px ui-monospace,monospace;color:#6c86a1}} .node-status{{position:absolute;right:12px;top:12px;width:10px;height:10px;border-radius:50%;background:var(--grey)}}
.node.completed .node-status,.completed.status-badge{{background:var(--green)}} .node.running .node-status,.running.status-badge{{background:var(--cyan);box-shadow:0 0 12px var(--cyan)}} .node.failed .node-status,.failed.status-badge{{background:var(--red)}} .node.pending{{opacity:.65}} .family-panel{{display:none}} .family-panel.active{{display:block}} .family-heading{{display:flex;justify-content:space-between;align-items:flex-start}} .status-badge{{border-radius:99px;padding:5px 10px;text-transform:uppercase;font-size:10px;font-weight:800;color:#06121e;background:var(--grey)}}
.table-wrap{{overflow:auto;border:1px solid var(--line);border-radius:10px}} table{{border-collapse:collapse;width:100%;min-width:700px}} th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid #20364f}} th{{color:#8fa9c3;background:#0b1828;font-size:11px;text-transform:uppercase;letter-spacing:.06em}} td{{font-family:ui-monospace,monospace;font-size:12px}} .pill{{padding:3px 7px;border-radius:99px;background:#31445a;font-size:10px}} .pill.completed,.pill.cached{{background:#123f37;color:#7ff0c1}} .pill.failed,.pill.aborted{{background:#4d1f2b;color:#ff93a4}} figure{{margin:18px 0;background:#fff;border-radius:12px;padding:10px}} figure img{{display:block;max-width:100%;margin:auto}} figcaption{{color:#50647b;padding:8px 4px 2px}} .empty{{color:var(--muted);border:1px dashed var(--line);border-radius:10px;padding:20px}} footer{{color:#607994;font-size:11px;margin:20px 2px}}
.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}} .gallery figure{{margin:0;min-width:0}} .failure-evidence{{background:#24131c;border:1px solid #733044;border-radius:16px;padding:18px;margin-top:14px;box-shadow:0 18px 55px #0004}} details{{background:#120f18;border:1px solid #4c2936;border-radius:10px;margin-top:10px;padding:10px 12px}} summary{{cursor:pointer;font-weight:700;color:#ff9bab}} pre{{white-space:pre-wrap;word-break:break-word;max-height:340px;overflow:auto;background:#080d16;border-radius:8px;padding:12px;color:#d8e5f5;font:11px/1.45 ui-monospace,monospace}} .evidence-note{{color:var(--amber)}}
@media(max-width:700px){{.shell{{padding:15px}}header{{display:block}}.live{{margin-top:12px;width:max-content}}}}
</style></head><body><div class="shell">
<header><div><span class="eyebrow">CRISPR Pipeline · execution dashboard</span><h1>{html.escape(args.run_name)}</h1><div class="run-id">{html.escape(args.run_id)}</div></div><div class="live {html.escape(args.status.lower())}"><i></i><span>{html.escape(args.status.upper())}</span></div></header>
<div class="summary">{metric_card("Completed", counts["COMPLETED"] + counts["CACHED"], "tasks")}{metric_card("Failed", counts["FAILED"] + counts["ABORTED"], "tasks")}{metric_card("Cached", counts["CACHED"], "tasks")}{metric_card("Task runtime", fmt_seconds(total_runtime), "aggregate")}{metric_card("Guides", guide.get("row_count", "—"), "validated")}</div>
<div class="graph-card"><div class="graph-head"><div><span class="eyebrow">Live dependency view</span><h2>Pipeline execution</h2></div><p>Click a family to inspect its QC and tasks</p></div><nav class="graph">{"".join(graph_nodes)}</nav></div>
{failures}{"".join(sections)}<footer>Generated {generated} · Self-contained W&amp;B HTML media · No credentials, FASTQs or unbounded task logs embedded</footer></div>
<script>function selectFamily(id){{document.querySelectorAll('.node,.family-panel').forEach(x=>x.classList.remove('active'));document.querySelector('[data-family="'+id+'"]').classList.add('active');document.getElementById('family-'+id).classList.add('active');}}selectFamily('{current}');</script>
</body></html>'''


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--status", default="running")
    parser.add_argument("--guide-report", type=Path, required=True)
    parser.add_argument("--seqspec-table", type=Path, required=True)
    parser.add_argument("--seqspec-image", type=Path)
    parser.add_argument("--qc-metrics-json", type=Path)
    parser.add_argument("--artifact-dir", type=Path)
    parser.add_argument("--nextflow-log", type=Path)
    parser.add_argument("--tail-lines", type=int, default=30)
    parser.add_argument("--max-image-bytes", type=int, default=10_000_000)
    parser.add_argument("--max-html-bytes", type=int, default=20_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    document = render(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document, encoding="utf-8")
    if args.output.stat().st_size > args.max_html_bytes:
        args.output.unlink()
        raise SystemExit(f"dashboard exceeds --max-html-bytes ({args.max_html_bytes})")
    print(args.output)
    print(f"dashboard_bytes={args.output.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
