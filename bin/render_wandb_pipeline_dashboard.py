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
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dashboard_presentation import enhance_dashboard_html


FAMILIES = [
    ("input", "Input validation and preparation", "Samples, metadata, references and covariates"),
    ("seqspec", "Library structure QC", "SeqSpec coordinates and guide capture"),
    ("mapping", "Library mapping", "Gene expression (GEX), guides and optional hashtags"),
    ("preprocessing", "GEX preprocessing", "Cell calling and per-set RNA quality filters"),
    ("mudata", "MuData assembly", "Align barcodes across modalities"),
    ("guide_assignment", "Guide assignment", "Guide-to-cell calls"),
    ("postconcat_qc", "Post-concatenation QC", "MT, clone filtering, normalization, PCA and UMAP"),
    ("inference", "Inference", "Estimate effects on gene expression"),
    ("evaluation", "Evaluation", "Controls and benchmarking"),
    ("final", "Final dashboard", "Published report and artifacts"),
]

CHECKPOINT_FAMILIES = {
    "input", "seqspec", "mapping", "preprocessing", "mudata", "guide_assignment", "postconcat_qc",
}

CATEGORY_FLOWS = {
    "input": [
        ("Validate inputs", "Samplesheet paths, modalities, measurement sets, guide metadata and provenance"),
        ("Prepare references", "Download or reuse genome/GTF resources and build guide/hash references"),
        ("Resolve covariates", "Prepare the covariates and batch fields requested for inference"),
    ],
    "seqspec": [
        ("Read supplied SeqSpec", "Use the declared barcode, UMI and library regions to configure mapping"),
        ("Check sampled reads", "Scan read mates and orientations for guide matches to check where the guide signal is found"),
        ("Report the evidence", "Show the best observed configuration per read pair and the plots used to assess it"),
    ],
    "mapping": [
        ("Quantify GEX", "Use kb/kallisto-bustools with the configured chemistry to create an unfiltered RNA count matrix"),
        ("Quantify feature barcodes", "Map guide reads and optional hashtag reads against their references to create count matrices"),
        ("Assemble matrices", "Combine count AnnData (.h5ad) outputs within each modality"),
        ("Check sequencing saturation", "When enabled, use RNA mapping records and RNA-QC-retained barcodes to estimate GEX saturation"),
    ],
    "preprocessing": [
        ("Call cells per measurement set", "Apply the configured RNA barcode caller to each measurement set"),
        ("Apply RNA cell filters", "Apply the configured UMI floor, two-sided RNA UMI and detected-gene bounds, and optional Scrublet"),
        ("Concatenate retained cells", "Combine independently filtered measurement-set matrices"),
        ("Pass retained RNA cells", "Carry raw RNA counts into modality intersection and guide assignment"),
    ],
    "mudata": [
        ("Intersect barcodes", "Match RNA-QC-retained cells to mapped guide barcodes and include the hash modality when enabled"),
        ("Assemble modalities", "Create the shared MuData object and preserve QC/provenance fields"),
    ],
    "guide_assignment": [
        ("Prepare assignment", "Select the configured capture and assignment model"),
        ("Call guide-positive cells", "Convert guide UMI evidence into binary cell–guide assignments with the configured model"),
        ("Concatenate and filter", "Combine per-set calls, then apply the configured assigned-guide eligibility and count rules"),
    ],
    "postconcat_qc": [
        ("Qualified cell intersection", "Start with GEX-qualified cells, assigned guides, and filtered HTO singlets when enabled"),
        ("Filter mitochondrial counts", "Apply the mitochondrial cutoff before temporary normalization, PCA and UMAP"),
        ("Recompute if clones are removed", "Repeat mitochondrial filtering and embedding QC on the retained raw counts"),
        ("Finalize MuData", "Apply gene-support filtering and deliver raw counts to inference"),
    ],
    "inference": [
        ("Define tests", "Build intended/local guide–gene pairs and global tests from validated metadata"),
        ("Fit methods", "Run configured SCEPTRE and/or PerTurbo local and global analyses"),
        ("Check calibration", "Compare SCEPTRE and PerTurbo against non-targeting pseudo-elements when enabled"),
        ("Assemble results", "Combine result chunks and write separate guide–gene and element–gene catalogs"),
    ],
    "evaluation": [
        ("Evaluate controls", "Compare intended effects, non-targeting controls and benchmark truth sets"),
    ],
    "final": [
        ("Collect outputs", "Gather available QC, inference, evaluation and provenance from stages 1–9"),
        ("Build report", "Create the complete local pipeline dashboard"),
        ("Publish the report", "Save the local HTML and refresh Weights & Biases when enabled"),
    ],
}


def family_for(process: str) -> str:
    value = process.lower()
    leaf = value.split(":")[-1]
    if leaf in {"save_checkpoint_dashboard", "stage_checkpoint_dashboard", "render_continuation_dashboard", "preinference_checkpoint"}:
        return "final"
    if "ntc_calibration" in value or "resolve_inference_mode" in value:
        return "inference"
    if any(key in leaf for key in ('embedding_before_clone', 'embedding_after_clone',
                                   'postconcat_embedding_qc', 'remove_clonal_cells', 'filter_hto_post_clone')):
        return 'postconcat_qc'
    if "seqspec" in value:
        return "seqspec"
    if "guide_mapping_qc" in leaf:
        return "mapping"
    if "sequencing_saturation" in leaf:
        return "mapping"
    if any(key in leaf for key in (
        "downloadreference", "skipgenomedownload", "skipgtfdownload",
        "createguideref", "createhashingref", "prepare_covariate",
    )):
        return "input"
    if any(key in value for key in ("mapping_rna_pipeline", "mapping_guide_pipeline", "mapping_hashing_pipeline")):
        if any(key in leaf for key in ("downloadreference", "seqspecparser", "createguideref", "createhashingref")):
            return "input"
        return "mapping"
    if any(key in value for key in ("preprocessing_pipeline", "preprocessanndata", "doublets", "filter_hashing")):
        return "preprocessing"
    if "guide_assignment" in value or "prepare_assignment" in value:
        return "guide_assignment"
    if any(key in leaf for key in ("createmudata", "anndata_concat", "hashing_concat")):
        return "mudata"
    if any(key in value for key in ("inference", "sceptre_chunk", "perturbo", "mergedresults", "catalog", "mergemudata")):
        return "inference"
    if any(key in value for key in (
        "evaluation", "additional_qc", "benchmark",
        "remove_clonal_cells", "clone_removal",
    )):
        return "evaluation"
    if "dashboard" in value or "publishfiles" in value:
        return "final"
    return "input"


def read_trace(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if not row.get("process"):
            row["process"] = (row.get("name") or "").split(" (", 1)[0]
    return rows


TASK_HANDLER_RE = re.compile(
    r"TaskHandler\[id:\s*(?P<task_id>\d+);\s*name:\s*(?P<name>.*?);\s*"
    r"status:\s*(?P<status>[A-Z]+);\s*exit:\s*(?P<exit>.*?);\s*"
    r"error:\s*.*?;\s*workDir:\s*(?P<workdir>[^\]]+)\]"
)
SUBMITTED_PROCESS_RE = re.compile(r"(?m)Submitted process > (?P<name>.+)$")


def merge_live_tasks(rows: list[dict[str, str]], nextflow_log: Path | None) -> list[dict[str, str]]:
    """Add the latest still-running TaskHandlers omitted from Nextflow trace.tsv."""
    if not nextflow_log or not nextflow_log.exists():
        return rows
    text = nextflow_log.read_bytes()[-2_000_000:].decode("utf-8", "replace")
    latest: dict[str, dict[str, str]] = {}
    for match in TASK_HANDLER_RE.finditer(text):
        item = match.groupdict()
        item["process"] = item["name"]
        item["duration"] = "in progress"
        item["realtime"] = "in progress"
        item["peak_rss"] = "—"
        latest[item["task_id"]] = item
    recorded = {row.get("task_id", "") for row in rows}
    recorded_names = {row.get("name", "") for row in rows}
    live = [
        item for task_id, item in latest.items()
        if task_id not in recorded
        and item["name"] not in recorded_names
        and item["status"] in {"NEW", "SUBMITTED", "RUNNING"}
    ]
    live_names = {item["name"] for item in live}
    submitted = {}
    for match in SUBMITTED_PROCESS_RE.finditer(text):
        name = match.group("name").strip()
        submitted[name] = {
            "task_id": f"submitted:{name}",
            "name": name,
            "process": name,
            "status": "SUBMITTED",
            "exit": "-",
            "workdir": "",
            "duration": "in progress",
            "realtime": "in progress",
            "peak_rss": "—",
        }
    live.extend(
        item for name, item in submitted.items()
        if name not in recorded_names and name not in live_names
    )
    return rows + live


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
    elif statuses["RUNNING"] + statuses["SUBMITTED"] + statuses["NEW"]:
        status = "running"
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
        "running": statuses["RUNNING"] + statuses["SUBMITTED"] + statuses["NEW"],
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


def enabled(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "yes", "1"}


def total_counts(rows: list[dict[str, Any]], key: str) -> int | None:
    values = [number(row.get(key)) for row in rows]
    if values and all(value is not None and value >= 0 and value.is_integer() for value in values):
        return sum(int(value) for value in values)
    return None


def modality_label(value: Any) -> Any:
    return {"trans": "GEX", "scrna": "GEX", "rna": "GEX", "grna": "Guide",
            "guide": "Guide", "hash": "HTO", "hashing": "HTO"}.get(str(value).lower(), value)


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


def searchable_table(
    rows: list[dict[str, Any]], columns: list[tuple[str, str]], table_id: str, limit: int = 25
) -> str:
    if not rows:
        return '<div class="empty">Inference results are not available yet.</div>'
    table = data_table(rows, columns, limit=limit).replace("<table>", f'<table id="{html.escape(table_id)}">', 1)
    return (
        f'<input class="table-search" type="search" placeholder="Search the {len(rows[:limit])} displayed rows" '
        f'aria-label="Search displayed guide–gene pairs" aria-controls="{html.escape(table_id)}">' + table
    )


def overall_row(section: dict[str, Any]) -> dict[str, Any]:
    rows = section.get("rows", []) if isinstance(section, dict) else []
    return next((row for row in rows if row.get("batch") == "all"), rows[0] if rows else {})


def parameter_stage(name: str) -> str | None:
    if name.startswith(("CLONE_REMOVAL_", "QC_EMBEDDING_")) or name in {
        "QC_pct_mito", "QC_min_cells_per_gene", "ENABLE_CLONE_REMOVAL", "ENABLE_POSTCONCAT_EMBEDDING_QC",
    }:
        return "Post-concatenation QC"
    if name.startswith("GUIDE_ASSIGNMENT_") or name == "QC_require_assigned_guide":
        return "Guide assignment"
    if name.startswith("INFERENCE_"):
        return "Inference"
    if name.startswith("SATURATION_") or name == "ENABLE_SEQUENCING_SATURATION":
        return "Mapping"
    if name == "ENABLE_DATA_HASHING":
        return "Mapping, MuData assembly and post-concatenation QC"
    if name == "QC_batch_col":
        return "RNA, guide and inference QC"
    if name == "Multiplicity_of_infection":
        return "MuData assembly and guide assignment"
    if name.startswith(("QC_", "SCRUBLET_")) or name == "ENABLE_SCRUBLET":
        return "GEX preprocessing"
    return None


def settings_content(params: dict[str, Any]) -> str:
    schema_path = Path(__file__).resolve().parent.parent / "nextflow_schema.json"
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        schema = {}
    definitions = [schema, *schema.get("$defs", {}).values()]
    properties = {key: prop for group in definitions for key, prop in group.get("properties", {}).items()}
    required = {key for group in definitions for key in group.get("required", [])}
    rows = []
    for key, value in sorted(params.items()):
        prop = properties.get(key, {})
        kind = prop.get("type")
        rows.append({
            "parameter": key, "value": value,
            "type": " / ".join(kind) if isinstance(kind, list) else kind,
            "optional": key not in required if prop else None,
            "description": prop.get("description"), "stage": parameter_stage(key),
        })
    return "<h3>Resolved analysis settings</h3>" + data_table(rows, [
        ("parameter", "Parameter"), ("value", "Value"), ("type", "Type"),
        ("optional", "Optional"), ("description", "Description"), ("stage", "Pipeline stage"),
    ], limit=max(1, len(rows)))


def run_cohort_flow(data: dict[str, Any]) -> str:
    observed = data.get("observed_metrics", {})
    additional = observed.get("additional_qc", {})
    rna = observed.get("measurement_set_rna_qc", {}).get("rows", [])
    overlap = additional.get("guide_mapping_qc", {}).get("rows", [])
    assignment = overall_row(additional.get("guide", {}))
    assignment_filter = additional.get("guide_assignment_filter", {}).get("rows", [])
    total_filter = next((row for row in assignment_filter if row.get("measurement_set") == "all"), {})
    summary = data.get("dashboard_summary", {}).get("filtering_summary", {})
    final = summary.get("final_mudata", {})
    stages = [
        ("Input RNA barcodes", total_counts(rna, "input_barcodes")),
        ("After RNA cell calling", total_counts(rna, "post_knee_cells")),
        ("After RNA QC", total_counts(rna, "retained_cells")),
        ("Exact RNA–guide matches", total_counts(overlap, "overlap_cells")),
        ("Cells with an assigned guide", assignment.get("n_cells_with_guide")
         if assignment.get("assigned_guide_counts_source") == "mudata_concat" else None),
        ("After guide-count filter", total_filter.get("cells_after")),
        ("Final MuData cells", final.get("cells")),
    ]
    if not any(value is not None for _, value in stages):
        return ""
    steps = "".join(
        '<div class="flow-step"><strong>' + html.escape(label) + '</strong><span>'
        + html.escape(display_value(value)) + '</span></div>'
        for label, value in stages
    )
    guide_count = total_counts(overlap, "guide_cells")
    guide_note = (
        '<p>Guide mapping contributed ' + html.escape(display_value(guide_count))
        + ' mapped guide barcodes before intersection.</p>' if guide_count is not None else ""
    )
    hto_note = (
        '<p>HTO qualification applies after guide assignment. <a href="#family-postconcat_qc" '
        'data-open-family="postconcat_qc">View HTO and downstream QC</a>.</p>'
        if enabled(data.get("parameters", {}).get("selected_qc_and_inference_params", {}).get("ENABLE_DATA_HASHING")) else ""
    )
    return ('<h3>Cell counts across stages</h3><p>Counts come from the stated stage. '
            'A dash indicates a count unavailable in this report.</p>' + guide_note
            + '<div class="flow-steps">' + steps + '</div>' + hto_note)


def category_flow(data: dict[str, Any], family: str) -> str:
    steps = CATEGORY_FLOWS.get(family, [])
    if not steps:
        return ""
    nodes = "".join(
        '<div class="flow-step"><strong>' + html.escape(title) + '</strong>'
        '<span>' + html.escape(description) + '</span></div>'
        for title, description in steps
    )
    resolved = ""
    params = data.get("parameters", {}).get("selected_qc_and_inference_params", {})
    if family == "preprocessing" and data:
        fields = [
            ("Barcode caller", params.get("QC_barcode_filter")),
            ("Minimum RNA UMI", params.get("QC_min_counts_per_cell")),
            ("RNA UMI MAD", params.get("QC_MAD_total_counts")),
            ("Detected-gene MAD", params.get("QC_MAD_n_genes")),
            ("Scrublet", params.get("ENABLE_SCRUBLET")),
        ]
        if enabled(params.get("ENABLE_SCRUBLET")):
            fields.extend([
                ("Scrublet profile", params.get("SCRUBLET_assay_type")),
                ("Scrublet rate override", params.get("SCRUBLET_expected_doublet_rate")),
                ("Scrublet PCA components", params.get("SCRUBLET_n_prin_comps")),
                ("Scrublet adaptive PCA fallback", params.get("SCRUBLET_adaptive_pca_fallback")),
            ])
        resolved = '<h4>Resolved filter values</h4><div class="filter-chips">' + "".join(
            '<span><b>' + html.escape(label) + ':</b> ' + html.escape(display_value(value)) + '</span>'
            for label, value in fields if value is not None
        ) + '</div><p class="flow-note">The RNA-count and detected-gene MAD filters are two-sided. '
        'The configured barcode caller, RNA UMI floor, MAD filters and optional Scrublet run per measurement set. '
        'Counts here stop at per-set RNA QC. Later filters and their counts appear in their own sections.</p>'
    if family == "postconcat_qc" and data:
        fields = [("Maximum mito %", params.get("QC_pct_mito")),
                  ("Minimum gene cell fraction", params.get("QC_min_cells_per_gene"))]
        resolved = '<h4>Downstream filter values</h4><div class="filter-chips">' + "".join(
            '<span><b>' + html.escape(label) + ':</b> ' + html.escape(display_value(value)) + '</span>'
            for label, value in fields if value is not None
        ) + '</div><p>These filters follow guide assignment; final-cohort statistics below include downstream cell QC.</p>'
    if family in {"guide_assignment", "inference"}:
        if family == "guide_assignment":
            method = str(params.get("GUIDE_ASSIGNMENT_method", "")).lower()
            keys = [("Assignment method", "GUIDE_ASSIGNMENT_method"),
                    ("Capture method", "GUIDE_ASSIGNMENT_capture_method")]
            if method == "sceptre":
                keys.extend([("Assignment probability", "GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold"),
                             ("EM repetitions", "GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep")])
            if method == "cleanser":
                keys.append(("Assignment probability", "GUIDE_ASSIGNMENT_cleanser_probability_threshold"))
        else:
            keys = [("Inference method", "INFERENCE_method"),
                    ("Target pairing", "INFERENCE_target_guide_pairing_strategy"),
                    ("Maximum local distance (bp)", "INFERENCE_max_target_distance_bp")]
        chips = "".join(
            '<span><b>' + label + ':</b> ' + html.escape(display_value(params[key])) + '</span>'
            for label, key in keys if params.get(key) is not None
        )
        if chips:
            resolved = '<h4>Resolved settings</h4><div class="filter-chips">' + chips + '</div>'
    return '<div class="process-flow"><h3>Processing steps</h3><div class="flow-steps">' + nodes + '</div>' + resolved + '</div>'


def qc_metrics_content(data: dict[str, Any], family: str, root: Path | None = None) -> str:
    if not data:
        return ""
    observed = data.get("observed_metrics", {})
    additional = observed.get("additional_qc", {})
    summary = data.get("dashboard_summary", {}).get("filtering_summary", {})

    if family == "input":
        params = data.get("parameters", {}).get("selected_qc_and_inference_params", {})
        return settings_content(params)

    if family == "mapping":
        rows = []
        for item in observed.get("mapping_json", []):
            metrics = item.get("metrics", {})
            rows.append({
                "measurement_set": item.get("measurement_set"), "modality": modality_label(item.get("modality")),
                "processed": metrics.get("n_processed"), "pseudoaligned": metrics.get("p_pseudoaligned"),
                "unique": metrics.get("p_unique"), "reads_onlist": metrics.get("percentageReadsOnOnlist"),
                "barcodes_onlist": metrics.get("percentageBarcodesOnOnlist"),
                "median_umis": metrics.get("medianUMIsPerBarcode"),
            })
        return ("<h3>Mapping and barcode capture</h3>"
                "<p>Each row describes one library. The on-list columns report the share of reads "
                "and observed barcodes matching the barcode whitelist. Median UMI/barcode comes from "
                "kb inspection of mapping records before RNA cell calling.</p>") + data_table(rows, [
            ("measurement_set", "Measurement set"), ("modality", "Modality"),
            ("processed", "Reads processed"), ("pseudoaligned", "Pseudoaligned %"),
            ("unique", "Unique %"), ("reads_onlist", "Reads on-list %"),
            ("barcodes_onlist", "Barcodes on-list %"), ("median_umis", "Median UMI/barcode"),
        ])

    if family == "preprocessing":
        rows = observed.get("measurement_set_rna_qc", {}).get("rows", [])
        cards = [
            metric_card("Input RNA barcodes", display_value(total_counts(rows, "input_barcodes")), "unfiltered RNA matrix rows"),
            metric_card("After barcode calling", display_value(total_counts(rows, "post_knee_cells")), "per-set total"),
            metric_card("After RNA UMI floor", display_value(total_counts(rows, "post_min_counts_cells")), "per-set total"),
            metric_card("After MAD policy", display_value(total_counts(rows, "post_mad_cells")), "per-set total"),
            metric_card("After RNA QC", display_value(total_counts(rows, "retained_cells")), "before modality intersection and guide assignment"),
        ]
        return '<h3>RNA preprocessing counts</h3><p>Input barcodes count rows of unfiltered RNA AnnData across measurement sets. Cell-calling and QC totals describe successive retained populations. <a href="#family-postconcat_qc" data-open-family="postconcat_qc">View pooled RNA, guide and available HTO QC</a>.</p><div class="metrics">' + "".join(cards) + "</div>" + data_table(
            rows, [("measurement_set", "Measurement set"), ("input_barcodes", "Input barcodes"),
            ("post_knee_cells", "After barcode calling"), ("post_min_counts_cells", "After UMI floor"),
            ("post_mad_cells", "After MAD policy"), ("retained_cells", "After RNA QC"),
            ("scrublet_status", "Scrublet status")], limit=max(1, len(rows)))

    if family == "mudata":
        rows = additional.get("guide_mapping_qc", {}).get("rows", [])
        return '<h3>RNA–guide modality intersection</h3><p>Counts at modality assembly, before assigned-guide and downstream cell filters.</p>' + data_table(
            rows, [("measurement_set", "Measurement set"), ("rna_cells", "RNA QC cells"),
                   ("guide_cells", "Guide-mapped cells"), ("overlap_cells", "Exact RNA–guide overlap")],
            limit=max(1, len(rows)))

    if family == "guide_assignment":
        overall = overall_row(additional.get("guide", {}))
        # Historical per-batch guide statistics describe final cells, not this stage.
        assigned = overall if overall.get("assigned_guide_counts_source") == "mudata_concat" else {}
        before = assigned.get("n_cells_before_assigned_guide_filter")
        after = assigned.get("n_cells_with_guide")
        rate = number(assigned.get("frac_cells_with_guide"))
        filters = additional.get("guide_assignment_filter", {}).get("rows", [])
        cards = [
            metric_card("Before assigned-guide filter", display_value(before), "assignment-stage cohort"),
            metric_card("Cells with an assigned guide", display_value(after), "before downstream cell QC"),
            metric_card("Cells without assigned guides", display_value(assigned.get("n_cells_without_assigned_guide")), "assignment-stage cohort"),
            metric_card("Assignment rate", f"{rate:.2%}" if rate is not None else "—", "with a guide / before assigned-guide filter"),
        ]
        content = ('<h3>Guide assignment and cell filtering</h3><p>Counts here precede mitochondrial filtering. '
                   '<a href="#family-postconcat_qc" data-open-family="postconcat_qc">View final-cohort guide statistics</a>.'
                   '</p><div class="metrics">' + "".join(cards) + '</div>')
        return content if not filters else content + data_table(
            filters, [("measurement_set", "Measurement set"), ("filter_label", "Assignment filter"),
                      ("cells_before", "Before"), ("cells_after", "After"), ("cells_removed", "Removed"),
                      ("threshold", "Threshold")], limit=max(1, len(filters)))

    if family == "postconcat_qc":
        gene = additional.get("gene", {})
        overall = overall_row(gene)
        final = summary.get("final_mudata", {})
        cards = [
            metric_card("Final cells", display_value(final.get("cells", overall.get("n_cells"))), "after downstream cell QC"),
            metric_card("Median RNA UMI", display_value(overall.get("umi_median")), "final cohort"),
            metric_card("Mean RNA UMI", display_value(overall.get("umi_mean")), "final cohort"),
            metric_card("Mean mito %", display_value(overall.get("mito_mean")), "final cohort"),
            metric_card("Gene features", display_value(final.get("gene_features")), "delivered MuData"),
            metric_card("Guide features", display_value(final.get("guide_features")), "delivered MuData"),
        ]
        content = '<h3>RNA quality after downstream cell QC</h3><div class="metrics">' + "".join(cards) + '</div>' + data_table(
            gene.get("rows", []), [("batch", "Batch"), ("n_cells", "Final cells"),
            ("umi_median", "Median UMI"), ("umi_mean", "Mean UMI"),
            ("mito_median", "Median mito %"), ("mito_mean", "Mean mito %")], limit=max(1, len(gene.get("rows", []))))
        guide = additional.get("guide", {})
        content += '<h3>Guide statistics in the final cohort</h3><p>These describe cells remaining after downstream QC, rather than the earlier assignment-stage cohort.</p>' + data_table(
            guide.get("rows", []), [("batch", "Batch"), ("n_cells", "Final cells"),
            ("guide_umi_median", "Median guide UMI"), ("n_cells_exactly_1_guide", "Exactly one guide"),
            ("guides_per_cell_mean", "Mean guides/cell")], limit=max(1, len(guide.get("rows", []))))
        return content

    if family == "evaluation":
        intended = overall_row(additional.get("intended_target", {}))
        global_qc = overall_row(additional.get("global_analysis", {}))
        tests = None
        summaries = sorted(root.rglob("global_analysis_per_guide_summary.tsv")) if root else []
        if summaries and all(path.read_bytes() == summaries[0].read_bytes() for path in summaries[1:]):
            with summaries[0].open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(rows) == number(global_qc.get("n_guides_tested")):
                tests = total_counts(rows, "n_tests")
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
        if tests is not None:
            cards.insert(5, metric_card("Valid trans tests", display_value(tests), "non-intended guide–gene pairs"))
        intended_columns = [
            ("n_guides_total", "Guides total"), ("n_guides_tested", "Guides tested"),
            ("fc_threshold", "FC threshold"), ("log2fc_threshold", "log2FC threshold"),
            ("pval_threshold", "P threshold"), ("n_strong_knockdowns", "Strong KD"),
            ("n_significant", "Significant"), ("n_strong_and_significant", "Strong + significant"),
            ("frac_strong_knockdowns", "Strong KD fraction"), ("frac_significant", "Significant fraction"),
            ("median_log2fc", "Median log2FC"), ("mean_log2fc", "Mean log2FC"),
            ("auroc", "AUROC"), ("auprc", "AUPRC"),
            ("n_eval_positives", "Eval positives"), ("n_eval_negatives", "Eval negatives"),
        ]
        global_columns = [
            ("n_guides_tested", "Guides tested"), ("n_targeting_guides", "Targeting"),
            ("n_non_targeting_guides", "Non-targeting"),
            ("median_significant_per_guide_targeting", "Median sig/targeting guide"),
            ("mean_significant_per_guide_targeting", "Mean sig/targeting guide"),
            ("median_significant_per_guide_nt", "Median sig/NT guide"),
            ("mean_significant_per_guide_nt", "Mean sig/NT guide"),
            ("total_significant_tests", "Significant tests"),
            ("median_genome_log2fc_targeting", "Median targeting log2FC"),
            ("median_genome_log2fc_nt", "Median NT log2FC"), ("fdr_method", "FDR method"),
            ("auroc", "AUROC"), ("auprc", "AUPRC"),
            ("n_validated_links", "Validated links"),
            ("n_eval_positives", "Eval positives"), ("n_eval_negatives", "Eval negatives"),
        ]
        return (
            '<h3>Effect and control evaluation</h3><div class="metrics">' + "".join(cards) + "</div>"
            '<h3>Intended-target QC details</h3><p>Counts guides with valid tests against their own intended gene. '
            'The table gives the fold-change and p-value rules used for this summary.</p>'
            + data_table([intended] if intended else [], intended_columns)
            + '<h3>Global-analysis QC details</h3><p>A trans test here checks a guide against a gene other '
              'than its own intended target. The totals include targeting and non-targeting guides with '
              'valid global results. The table identifies the multiple-testing correction.</p>'
            + data_table([global_qc] if global_qc else [], global_columns)
        )

    if family == "final":
        sources = data.get("sources", {})
        rows = [{"source": key, "value": value} for key, value in sorted(sources.items())]
        return (run_cohort_flow(data) + '<h3>Metric sources</h3><p>Inputs and recorded counts from the '
                f'run’s QC metrics file (schema {html.escape(str(data.get("schema_version", "—")))}).</p>'
                + data_table(rows, [("source", "Source"), ("value", "Value")]))
    return ""


def image_family(path: Path) -> str:
    # The output directory itself may contain "embedding_qc".  Classifying
    # the entire absolute path would therefore route every image to post-QC.
    value = path.as_posix().lower()
    parts = {part.lower() for part in path.parts}
    if "calibration" in parts or path.name.lower().startswith("ntc_"):
        return "inference"
    if any(term in parts for term in (
        'postconcat_embedding_qc', 'embedding_qc', 'embeddings', 'before_clone',
        'after_clone', 'clone_qc', 'clone_removal', 'clones', 'clone_filter',
        'hto_filter', 'hashing_qc',
    )):
        return 'postconcat_qc'
    if path.stem.lower() in {"cells_per_hto_barplot", "umap_hto", "umap_hto_singlets"}:
        return "postconcat_qc"
    if "additional_qc" in parts and ("gene" in parts or "guide" in parts):
        return "postconcat_qc"
    if "seqspec" in value:
        return "seqspec"
    if "guide_mapping_qc" in parts or "guide_mapping_orientation_qc" in value:
        return "mapping"
    if "sequencing_saturation" in value:
        return "mapping"
    if "benchmark_output" in parts or path.stem.lower() == "tf_benchmark":
        return "evaluation"
    if any(term in value for term in ("intended_target", "global_analysis", "evaluation", "volcano")):
        return "evaluation"
    if "figures" in parts and path.stem.lower() in {
        "gene_cells_per_batch", "gene_histograms", "gene_histograms_by_batch", "gene_knee_plot",
        "scrna_barcodes_umi_thresholds", "scrna_barcodes_detected_genes_thresholds",
        "guides_per_cell_histogram", "cells_per_guide_histogram", "guides_umi_thresholds", "guides_hist_num_sgrna",
    }:
        return "postconcat_qc"
    if any(term in value for term in ("guide_", "guides_", "sgrna", "cells_per_guide", "guides_per_cell", "hto_", "hashing_qc")):
        return "guide_assignment"
    if any(term in value for term in ("scrna", "rna_qc", "gene_", "knee_plot")):
        return "preprocessing"
    if any(term in value for term in ("loss_curve", "perturbo", "sceptre")):
        return "inference"
    return "final"


def embedding_stage(path: Path) -> str:
    for part in path.parts:
        if part in ('before_clone', 'after_clone'):
            return part
    metrics = path.parent / 'embedding_qc_metrics.json'
    if metrics.is_file():
        try:
            return str(json.loads(metrics.read_text()).get('stage', ''))
        except (OSError, ValueError):
            pass
    return ''


def collect_images(root: Path | None, max_bytes: int = 38_000_000,
                   omitted: list | None = None, checkpoint: Path | None = None) -> dict[str, list[tuple[Path, str]]]:
    selected: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    if not root or not root.exists():
        return selected
    seen: set[tuple[str, str]] = set()
    used = 0
    # Reserve the first slots for the central post-concatenation diagnostics.
    priorities = {'normalization_check.png', 'pca_variance_ratio.png', 'pca_qc_panel.png',
                  'pca_measurement_sets_colored.png', 'pca_by_measurement_set.png', 'umap_qc_panel.png',
                  'leiden_resolution_sweep.png', 'leiden_sweep_umap.png',
                  'cell_cycle_by_measurement_set.png', 'postconcat_qc_flow.png',
                  'ntc_calibration_qq_guide.png', 'ntc_calibration_qq_element.png',
                  'ntc_aucpr_curves_guide.png', 'ntc_aucpr_curves_element.png'}
    candidates = []
    for path in root.rglob("*.png"):
        relative = path.relative_to(root)
        if relative.parts[0] in {"preinference_checkpoint", "checkpoint_dashboard"}:
            continue
        family = image_family(relative)
        if checkpoint is None or family not in CHECKPOINT_FAMILIES:
            candidates.append((path, relative, family))
    if checkpoint is not None:
        for path in checkpoint.rglob("*.png"):
            relative = path.relative_to(checkpoint)
            candidates.append((path, relative, image_family(relative)))
    for path, relative, family in sorted(candidates, key=lambda item: (
        item[0].name not in priorities and "calibration" not in item[1].parts, str(item[0])
    )):
        size = path.stat().st_size
        content = path.read_bytes()
        # Preserve identical pre/post-clone images as distinct scientific views,
        # while collapsing duplicate published copies within the same stage.
        digest = (hashlib.sha256(content).hexdigest(), embedding_stage(path))
        if digest in seen:
            continue
        seen.add(digest)
        if used + size > max_bytes:
            if omitted is not None:
                omitted.append({'plot': str(relative), 'bytes': size,
                                'reason': 'image byte budget'})
            continue
        used += size
        encoded = base64.b64encode(content).decode("ascii")
        selected[family].append((path, encoded))
    return selected


def image_gallery(images: list[tuple[Path, str]], title: str = "QC visualizations") -> str:
    if not images:
        return ""
    figures = "".join(
        f'<figure><img loading="lazy" src="data:image/png;base64,{encoded}" alt="{html.escape((embedding_stage(path) + " " + path.stem).strip())}">'
        f'<figcaption>{html.escape((embedding_stage(path) + " · " + path.stem.replace("_", " ")).strip(" ·"))}</figcaption></figure>'
        for path, encoded in images
    )
    return f'<h3>{html.escape(title)}</h3><div class="gallery">{figures}</div>'


def measurement_set_names(data: dict[str, Any]) -> list[str]:
    observed = data.get("observed_metrics", {}) if data else {}
    names = {
        str(row.get("measurement_set"))
        for row in observed.get("mapping_json", [])
        if row.get("measurement_set")
    }
    names.update(
        str(row.get("measurement_set"))
        for row in observed.get("measurement_set_rna_qc", {}).get("rows", [])
        if row.get("measurement_set")
    )
    return sorted(names)


def _percent_bar(value: Any, label: str) -> str:
    numeric = number(value)
    width = max(0.0, min(100.0, numeric if numeric is not None else 0.0))
    shown = display_value(value)
    return (
        '<div class="mini-bar"><div class="mini-bar-label"><span>' + html.escape(label) +
        '</span><strong>' + html.escape(shown) + '%</strong></div>'
        f'<div class="mini-bar-track"><i style="width:{width:.2f}%"></i></div></div>'
    )


def mapping_measurement_set_cards(data: dict[str, Any]) -> str:
    observed = data.get("observed_metrics", {}) if data else {}
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in observed.get("mapping_json", []):
        if row.get("measurement_set"):
            grouped[str(row["measurement_set"])].append(row)
    if not grouped:
        return ""
    cards = []
    for measurement_set, rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: str(row.get("modality", "")))
        total_reads = sum(number(row.get("metrics", {}).get("n_processed")) or 0 for row in rows)
        modality_chips = "".join(
            '<span>' + html.escape(str(modality_label(row.get("modality", "unknown")))) + '</span>' for row in rows
        )
        modality_sections = []
        for row in rows:
            metrics = row.get("metrics", {})
            modality = str(modality_label(row.get("modality", "unknown")))
            facts = [
                ("Reads processed", metrics.get("n_processed")),
                ("Pseudoaligned", f'{display_value(metrics.get("p_pseudoaligned"))}%'),
                ("Unique", f'{display_value(metrics.get("p_unique"))}%'),
                ("Reads on-list", f'{display_value(metrics.get("percentageReadsOnOnlist"))}%'),
                ("Barcodes on-list", f'{display_value(metrics.get("percentageBarcodesOnOnlist"))}%'),
                ("Median UMI/barcode", metrics.get("medianUMIsPerBarcode")),
            ]
            fact_html = "".join(
                '<div class="measurement-fact"><span>' + html.escape(label) + '</span><strong>' +
                html.escape(display_value(value)) + '</strong></div>' for label, value in facts
            )
            bars = _percent_bar(metrics.get("p_pseudoaligned"), "Pseudoaligned")
            if metrics.get("percentageReadsOnOnlist") is not None:
                bars += _percent_bar(metrics.get("percentageReadsOnOnlist"), "Reads on-list")
            modality_sections.append(
                '<section class="modality-block"><h4>' + html.escape(modality) + '</h4>' +
                '<div class="measurement-facts">' + fact_html + '</div>' + bars + '</section>'
            )
        cards.append(
            '<details class="measurement-card"><summary><div><span class="measurement-name">' +
            html.escape(measurement_set) + '</span><span class="measurement-subtitle">' +
            html.escape(", ".join(str(modality_label(row.get("modality", "unknown"))) for row in rows)) +
            '</span></div><div class="measurement-summary"><strong>' +
            html.escape(display_value(int(total_reads))) + '</strong><span>processed reads across libraries</span></div>' +
            '<div class="modality-chips">' + modality_chips + '</div></summary>' +
            '<div class="measurement-body">' + "".join(modality_sections) + '</div></details>'
        )
    return (
        '<div class="measurement-section"><div class="measurement-section-head"><div>'
        '<span class="eyebrow">Measurement-set hierarchy</span><h3>Mapping by measurement set</h3>'
        '<p>Open a card for per-library rates. Each card total sums processed reads across libraries in that measurement set.</p></div>'
        f'<span class="measurement-count">{len(cards)} sets</span></div>' + "".join(cards) + '</div>'
    )


def guide_mapping_qc_content(root: Path | None) -> str:
    """Render configured orientation and post-mapping barcode recovery."""
    if not root or not root.exists():
        return ""
    reports = sorted(root.glob("**/guide_mapping_qc.json"))
    tables = sorted(root.glob("**/guide_mapping_qc.tsv"))
    if not reports:
        return ""
    try:
        report = json.loads(reports[-1].read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    overall = report.get("overall", {})
    orientation = report.get("configured_orientation", {})
    status = str(report.get("status", "UNKNOWN"))
    cards = [
        metric_card("Guide mapping recovery", status, "before guide assignment"),
        metric_card(
            "Guide orientation",
            "reverse complement" if orientation.get("reverse_complement_guides") else "as supplied",
            "guide reference",
        ),
        metric_card("Spacer tag", display_value(orientation.get("spacer_tag")), "guide search anchor"),
        metric_card("GEX cells", display_value(overall.get("rna_cells")), "after per-set RNA QC"),
        metric_card("Mapped guide barcodes", display_value(overall.get("guide_cells")), "before modality intersection"),
        metric_card("MuData barcodes", display_value(overall.get("mudata_intersection_cells")), "after modality intersection"),
        metric_card(
            "Recovered guide designs",
            f'{display_value(overall.get("expected_guides_with_nonzero_counts"))} / {display_value(overall.get("expected_guides"))}',
            "at least one mapped UMI",
        ),
    ]
    rows = []
    if tables:
        try:
            with tables[-1].open(newline="", encoding="utf-8", errors="replace") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        except OSError:
            rows = []
    columns = [
        ("measurement_set", "Measurement set"),
        ("rna_cells", "GEX cells after RNA QC"),
        ("guide_cells", "Mapped guide barcodes"),
        ("overlap_cells", "Exact overlap"),
        ("guide_to_rna_fraction", "Guide barcodes / GEX cells"),
        ("overlap_to_guide_fraction", "Overlap / guide barcodes"),
        ("overlap_to_rna_fraction", "Overlap / GEX cells"),
        ("status", "Status"),
        ("reason", "Reason"),
    ]
    failures = report.get("failures", [])
    failure_html = ""
    if failures:
        failure_html = '<div class="callout warn"><strong>Recovery warning</strong><br>' + html.escape(" | ".join(map(str, failures))) + "</div>"
    thresholds = report.get("thresholds", {})
    rules = [
        ("min_overlap_cells_per_set", "Exact overlaps per set"),
        ("min_guide_to_rna_fraction", "Guide barcodes / GEX cells"),
        ("min_overlap_to_rna_fraction", "Overlap / GEX cells"),
        ("min_overlap_to_guide_fraction", "Overlap / guide barcodes"),
        ("min_recovered_guide_fraction", "Expected designs recovered across sets"),
    ]
    published_rules = [{"rule": label, "minimum": thresholds[key]} for key, label in rules if key in thresholds]
    gate_details = ""
    if published_rules:
        gate_details = ('<details class="supporting-details"><summary>Recovery criteria used in this run</summary>'
                        + data_table(published_rules, [("rule", "Required measure"), ("minimum", "Minimum")])
                        + "</details>")
    return (
        '<div class="measurement-section"><div class="measurement-section-head"><div>'
        '<span class="eyebrow">Pre-assignment validation</span><h3>Guide mapping and orientation QC</h3>'
        '<p>Shows the orientation used for mapping and the recorded recovery check. The table compares '
        'RNA-QC-retained GEX cells with mapped guide barcodes before assignment. '
        'The guide-to-GEX ratio can exceed 1. Exact overlap uses identical barcode IDs.</p>'
        '</div></div><div class="metrics">' + "".join(cards) + "</div>" + failure_html +
        data_table(rows, columns, limit=max(1, len(rows))) + gate_details + "</div>"
    )


def preprocessing_measurement_set_cards(
    data: dict[str, Any], images: list[tuple[Path, str]]
) -> tuple[str, set[Path]]:
    observed = data.get("observed_metrics", {}) if data else {}
    rows = observed.get("measurement_set_rna_qc", {}).get("rows", [])
    if not rows:
        return "", set()
    used: set[Path] = set()
    cards = []
    plot_priority = ("knee_plot", "rna_qc_filter_flow", "rna_qc_filter_steps", "qc_distributions", "scrublet_scores")
    for row in sorted(rows, key=lambda item: str(item.get("measurement_set", ""))):
        measurement_set = str(row.get("measurement_set", "unknown"))
        matching = [(path, encoded) for path, encoded in images if measurement_set.lower() in path.name.lower()]
        matching.sort(key=lambda item: next(
            (index for index, prefix in enumerate(plot_priority) if prefix in item[0].name.lower()),
            len(plot_priority),
        ))
        used.update(path for path, _ in matching)
        input_cells = number(row.get("input_barcodes"))
        retained = number(row.get("retained_cells"))
        retained_pct = f"{100.0 * retained / input_cells:.2f}%" if input_cells and retained is not None else "—"
        scrublet_status = str(row.get("scrublet_status", "")).lower()
        show_scrublet = scrublet_status == "applied" or scrublet_status.startswith("skipped") or (
            not scrublet_status and enabled(row.get("scrublet_enabled"))
        )
        scrublet_label = "After Scrublet"
        if scrublet_status.startswith("skipped"):
            scrublet_label = "Scrublet skipped (fit error)" if scrublet_status == "skipped_error" else "Scrublet skipped"
        if not show_scrublet:
            matching = [(path, encoded) for path, encoded in matching if "scrublet_scores" not in path.name.lower()]
        caller_labels = {"none": "After barcode policy", "knee": "After knee", "knee2": "After knee2"}
        stages = [
            ("Input barcodes", row.get("input_barcodes")),
            (caller_labels.get(row.get("barcode_filter"), "After barcode calling"), row.get("post_knee_cells")),
            ("After RNA UMI floor", row.get("post_min_counts_cells")),
            ("After RNA UMI and gene MAD filters", row.get("post_mad_cells")),
        ]
        if show_scrublet:
            stages.append((scrublet_label, row.get("retained_cells")))
        stage_html = "".join(
            '<div class="filter-stage"><span>' + html.escape(label) + '</span><strong>' +
            html.escape(display_value(value)) + '</strong></div>' for label, value in stages
        )
        plots = "".join(
            f'<figure><img loading="lazy" src="data:image/png;base64,{encoded}" alt="{html.escape(path.stem)}">'
            f'<figcaption>{html.escape(path.stem.replace("_", " "))}</figcaption></figure>'
            for path, encoded in matching
        )
        facts = [] if not show_scrublet or scrublet_status.startswith("skipped") else [
            ("Scrublet removed", row.get("removed_by_scrublet")),
            ("Scrublet expected doublet rate", row.get("scrublet_expected_doublet_rate")),
            ("Scrublet PCA components used", row.get("scrublet_n_prin_comps_used")),
        ]
        fact_html = "".join(
            '<div class="measurement-fact"><span>' + html.escape(label) + '</span><strong>' +
            html.escape(display_value(value)) + '</strong></div>' for label, value in facts if value is not None
        )
        cards.append(
            '<details class="measurement-card preprocessing-card"><summary><div><span class="measurement-name">' +
            html.escape(measurement_set) + '</span></div>'
            '<div class="measurement-summary"><strong>' + html.escape(display_value(row.get("retained_cells"))) +
            f'</strong><span>RNA-QC retained · {retained_pct}</span></div></summary>'
            '<div class="measurement-body"><h4>Cell-filter sequence</h4><div class="filter-stages">' +
            stage_html + '</div><div class="measurement-facts">' + fact_html + '</div>' +
            ('<div class="measurement-plots">' + plots + '</div>' if plots else '<p class="empty">Plots have not been published for this set yet.</p>') +
            '</div></details>'
        )
    return (
        '<div class="measurement-section"><div class="measurement-section-head"><div>'
        '<span class="eyebrow">Measurement-set hierarchy</span><h3>Preprocessing by measurement set</h3>'
        '<p>Open a card to inspect per-set RNA filters and their retained cells. Counts precede guide assignment and downstream cell QC.</p></div>'
        f'<span class="measurement-count">{len(cards)} sets</span></div>' + "".join(cards) + '</div>', used
    )


def figure_measurement_set_cards(
    images: list[tuple[Path, str]],
    measurement_sets: list[str],
    title: str,
    description: str,
    card_subtitle: str,
) -> tuple[str, set[Path]]:
    cards = []
    used: set[Path] = set()
    for measurement_set in measurement_sets:
        matching = [
            (path, encoded)
            for path, encoded in images
            if measurement_set.lower() in path.name.lower()
        ]
        if not matching:
            continue
        used.update(path for path, _ in matching)
        plots = "".join(
            f'<figure><img loading="lazy" src="data:image/png;base64,{encoded}" alt="{html.escape(path.stem)}">'
            f'<figcaption>{html.escape(path.stem.replace("_", " "))}</figcaption></figure>'
            for path, encoded in matching
        )
        cards.append(
            '<details class="measurement-card figure-card"><summary><div><span class="measurement-name">' +
            html.escape(measurement_set) + '</span><span class="measurement-subtitle">' +
            html.escape(card_subtitle) + '</span></div><div class="measurement-summary"><strong>' +
            str(len(matching)) + '</strong><span>QC plot' + ('s' if len(matching) != 1 else '') +
            '</span></div></summary><div class="measurement-body"><div class="measurement-plots">' +
            plots + '</div></div></details>'
        )
    if not cards:
        return "", set()
    return (
        '<div class="measurement-section"><div class="measurement-section-head"><div>'
        '<span class="eyebrow">Measurement-set hierarchy</span><h3>' + html.escape(title) + '</h3><p>' +
        html.escape(description) + '</p></div><span class="measurement-count">' +
        str(len(cards)) + ' sets</span></div>' + "".join(cards) + '</div>',
        used,
    )


def measurement_set_filter_flow_content(root: Path | None) -> str:
    if not root or not root.exists():
        return ""
    candidates = sorted(root.glob("**/measurement_set_qc_filter_flow.tsv"))
    if not candidates:
        return ""
    with candidates[-1].open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [
        row for row in rows
        if not (row.get("filter_label") == "Scrublet doublet removal" and row.get("threshold") == "disabled")
    ]
    for row in rows:
        threshold = row.get("threshold", "")
        if threshold.startswith("rank ≤") and "RNA UMIs ≥" in threshold:
            row["threshold"] = threshold.replace("rank ≤ ", "Knee rank ", 1)
        for key in ("removed_percent", "retained_percent_of_input"):
            value = number(row.get(key))
            if value is not None:
                row[key] = f"{value:.2f}%"
    columns = [
        ("measurement_set", "Measurement set"),
        ("step_order", "Order"),
        ("filter_label", "Filter parameter"),
        ("threshold", "Resolved rule"),
        ("applied", "Applied"),
        ("cells_before", "Cells before"),
        ("cells_after", "Cells after"),
        ("cells_removed", "Removed"),
        ("removed_percent", "Removed %"),
        ("retained_percent_of_input", "Input retained %"),
    ]
    return (
        '<h3>Sequential per-measurement-set RNA filters</h3>'
        '<p>Removed % uses the input to that filter. Input retained % uses the original measurement-set '
        'barcode count. A UMI-driven knee call includes barcodes tied at the cutoff. An applied filter can '
        'remove zero cells when all pass.</p>'
        + data_table(rows, columns, limit=max(1, len(rows)))
    )


def assignment_filter_flow_content(root: Path | None, family: str) -> str:
    """Show the guide or HTO filters recorded for this stage."""
    if not root or not root.exists():
        return ""
    specifications = {
        "guide_assignment": (
            "**/guide_assignment_filter_flow.tsv",
            "Guide-assignment cell filter",
            "These limits are applied to the number of assigned guides per cell after guide calling. Counts removed are cells.",
        ),
        "postconcat_qc": (
            "**/hto_filter_flow.tsv",
            "HTO intersection filters",
            "In the embedding workflow HTO support and singlet retention are calculated after guide QC, before the parallel clone/embedding branches.",
        ),
    }
    columns = [
        ("measurement_set", "Measurement set"), ("step_order", "Order"),
        ("filter_label", "Filter parameter"), ("threshold", "Resolved rule"),
        ("applied", "Applied"), ("cells_before", "Cells before"),
        ("cells_after", "Cells after"), ("cells_removed", "Removed"),
        ("removed_percent", "Removed %"),
        ("retained_percent_of_input", "Input retained %"),
    ]
    pattern, title, note = specifications[family]
    candidates = sorted(root.glob(pattern))
    if not candidates:
        return ""
    rows = []
    seen_tables = set()
    for candidate in candidates:
        digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if digest in seen_tables:
            continue
        seen_tables.add(digest)
        with candidate.open(newline="", encoding="utf-8", errors="replace") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return f"<h3>{html.escape(title)}</h3><p>{html.escape(note)}</p>" + data_table(rows, columns, limit=max(1, len(rows)))


def postconcat_qc_content(root: Path | None) -> str:
    if not root or not root.exists():
        return '<p>Post-concatenation QC artifacts have not been published yet.</p>'
    stages = {}
    for path in sorted(root.rglob('embedding_qc_metrics.json')):
        try:
            row = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        stages[row.get('stage', str(path))] = (row, path.parent)
    if not stages:
        clone_paths = sorted(root.rglob('clone_metrics.tsv'))
        if clone_paths:
            with clone_paths[-1].open() as handle:
                clone_rows = list(csv.DictReader(handle, delimiter='\t'))
            if clone_rows:
                content = '<h3>Clone calling/removal</h3>' + data_table(
                    clone_rows, [(key, key.replace('_', ' ')) for key in clone_rows[0]]
                )
                content += '<p>Embedding QC was not produced by this run.</p>'
                return content
        return '<p>Post-concatenation QC artifacts have not been published yet.</p>'
    rows = []
    for _, (source, _) in sorted(stages.items(), reverse=True):
        row = source.copy()
        cell_cycle = source.get('cell_cycle', {})
        leiden = source.get('leiden', {})
        row.update(cell_cycle_status=cell_cycle.get('status', '—'),
                   cell_cycle_overlap=(f"S {cell_cycle.get('s_markers_present', '—')} / "
                                       f"G2M {cell_cycle.get('g2m_markers_present', '—')}"),
                   leiden_status=leiden.get('status', '—'),
                   leiden_diagnostic_clusters=leiden.get('diagnostic_clusters', '—'))
        rows.append(row)
    content = '<h3>Actual normalization and embedding settings</h3>' + data_table(rows, [
        ('stage', 'Stage'), ('status', 'Status'), ('input_cells', 'Input cells'),
        ('retained_cells', 'After MT'), ('cells', 'Cells embedded'), ('mito_cutoff', 'MT cutoff %'),
        ('normalization_target_sum', 'Median-depth target'), ('feature_selection', 'Feature selection'), ('genes', 'Expressed genes'),
        ('pca_features', 'PCA features'), ('requested_pcs', 'Requested PCs'), ('effective_pcs', 'Used PCs'),
        ('effective_neighbors', 'Neighbors'), ('genes_after', 'Delivered genes'),
        ('cell_cycle_status', 'Cell cycle'), ('cell_cycle_overlap', 'CC markers represented'),
        ('leiden_status', 'Leiden sweep'), ('leiden_diagnostic_clusters', 'Diagnostic clusters'),
        ('normalized_matrix_saved', 'Normalized matrix saved'), ('omitted_covariates', 'Unavailable covariates')])
    content += '<p>No batch correction. Independent pre/post-clone UMAPs can rotate; compare covariate patterns, not absolute positions. TAP-seq measures a targeted panel, not the entire transcriptome.</p>'
    for stage, (row, directory) in sorted(stages.items(), reverse=True):
        flow = directory / 'measurement_filter_flow.tsv'
        if flow.exists():
            with flow.open() as handle:
                flow_rows = list(csv.DictReader(handle, delimiter='\t'))
            content += f'<h3>{html.escape(stage)}: mitochondrial filtering by measurement set</h3>' + data_table(flow_rows, [
                ('measurement_set', 'Measurement set'), ('cells_before', 'Before'), ('cells_after', 'After'),
                ('cells_removed', 'Removed'), ('threshold', 'MT cutoff %')])
    clone_paths = sorted(root.rglob('clone_metrics.tsv'))
    if clone_paths:
        with clone_paths[-1].open() as handle:
            clone_rows = list(csv.DictReader(handle, delimiter='\t'))
        if clone_rows:
            content += '<h3>Clone calling/removal</h3>' + data_table(clone_rows, [(k, k.replace('_', ' ')) for k in clone_rows[0]])
            for row in clone_rows:
                if row.get('applicability') == 'low_power_warning':
                    content += '<p class="empty"><strong>Clone-calling applicability warning:</strong> ' + html.escape(row.get('applicability_note', '')) + '</p>'
    content += ('<h3>Scope relative to the reference QC scripts</h3>'
                '<p>Implemented: median-depth normalization, log1p, HVG/all-panel feature selection, '
                'scale clipping at 10, PCA variance, PCA/UMAP covariate panels, PCA measurement-set facets, '
                'cell-cycle scoring for eligible whole-transcriptome assays, native-igraph Leiden sweeps and '
                'measurement-set composition views. Cell-cycle scoring is skipped automatically for targeted '
                'TAP-seq panels unless explicitly forced. These are descriptive QC views; no transformed matrix, '
                'embedding or cluster assignment is saved into the inference MuData.</p>')
    return content


def clean_html_cell(value: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", value)).strip()


def final_inference_content(path: Path | None, limit: int = 25) -> str:
    """Mirror bounded top inference rows already selected by the final dashboard."""
    if not path or not path.is_file():
        return '<h3>Inference results</h3><div class="empty">Result tables will appear after the final dashboard selects the lowest-p-value pairs.</div>'
    source = path.read_text(encoding="utf-8", errors="replace")
    marker = re.compile(r'<h3[^>]*>Guide Inference</h3>', re.I)
    sections = []
    preferred = [
        "gene_name", "gene_id", "guide_id", "intended_target_name",
        "genomic_element_id", "genomic_element_name", "genomic_element_type",
        "genomic_element_chr", "genomic_element_start", "genomic_element_end",
        "guide_ids", "num_guides", "nPerturbedCells",
        "sceptre_log2_fc", "sceptre_p_value", "sceptre_q_value",
        "perturbo_log2_fc", "perturbo_p_value", "perturbo_q_value",
    ]
    for index, match in enumerate(marker.finditer(source)):
        start = match.end()
        table_start = source.find("<table", start)
        if table_start < 0:
            continue
        opening_end = source.find(">", table_start) + 1
        table_end = source.find("</table>", opening_end)
        if table_end < 0:
            continue
        cursor = opening_end
        for _ in range(limit + 1):
            found = source.find("</tr>", cursor, table_end)
            if found < 0:
                break
            cursor = found + len("</tr>")
        fragment = source[opening_end:cursor]
        columns = [clean_html_cell(cell) for cell in re.findall(r"<th[^>]*>(.*?)</th>", fragment, re.I | re.S)]
        rows = []
        for row_html in re.findall(r"<tr[^>]*>(.*?)</tr>", fragment, re.I | re.S):
            cells = [clean_html_cell(cell) for cell in re.findall(r"<td[^>]*>(.*?)</td>", row_html, re.I | re.S)]
            if cells and len(cells) == len(columns):
                rows.append(dict(zip(columns, cells)))
        label_match = re.search(r"<p[^>]*>\s*(Local Analysis|Global Analysis)\s*</p>", source[start:table_start], re.I)
        label = label_match.group(1) if label_match else f"Inference {index + 1}"
        labels = {"intended_target_name": "Intended target (as recorded)", "nPerturbedCells": "Perturbed cells"}
        selected = [(column, labels.get(column, column.replace("_", " "))) for column in preferred if column in columns]
        sections.append(
            f'<div class="result-block"><h3>{html.escape(label)}: top guide–gene pairs</h3>'
            f'<p>Preview of {len(rows)} lowest-p-value rows selected by the final report. Search covers these displayed rows.</p>'
            + searchable_table(rows, selected, f"inference-{index}", limit=limit) + "</div>"
        )
    if not sections:
        return '<h3>Inference results</h3><div class="empty">No Guide Inference tables were found in the final dashboard.</div>'
    return "".join(sections)


def calibration_artifact_content(root: Path | None) -> str:
    """Show the PR #103 NTC check beside the inference methods it calibrates."""
    if not root:
        return ""
    candidates = (
        root / "pipeline_outputs" / "calibration",
        root / "pipeline_dashboard" / "benchmark_output" / "calibration",
    )
    calibration_dir = next((candidate for candidate in candidates if candidate.is_dir()), None)
    if calibration_dir is None:
        return ""
    heading = '<div class="result-block"><h3>Non-targeting calibration</h3>'
    summary_file = calibration_dir / "ntc_calibration_summary.json"
    skipped_file = calibration_dir / "ntc_calibration_skipped.txt"
    if skipped_file.is_file():
        return (heading + '<p>The calibration pass was skipped.</p><div class="qc-callout">'
                + html.escape(sanitized_tail(skipped_file, 12)) + '</div></div>')
    if not summary_file.is_file():
        return heading + '<p>The calibration report is still being prepared.</p></div>'
    try:
        summary = json.loads(summary_file.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return heading + '<p>The calibration summary could not be read yet.</p></div>'
    if summary.get("status") != "ok":
        return heading + '<p>The calibration check was not completed.</p></div>'
    calibration_rows = summary.get("calibration", [])
    auc_rows = summary.get("auprc", [])
    calibration_columns = [
        ("level", "Level"), ("pair_set", "Pair set"), ("method", "Method"),
        ("n_tested", "Tested"), ("frac_p_lt_0.05", "Fraction p < 0.05"),
        ("ks_statistic", "KS distance"),
    ]
    auc_columns = [
        (key, key.replace("_", " ").title())
        for key in ("level", "method", "negatives", "n_direct_target", "n_ntc_used", "positive_rate", "auprc", "auroc")
        if any(key in row for row in auc_rows)
    ]
    notes = summary.get("notes", [])
    note_text = "".join(f"<li>{html.escape(str(note))}</li>" for note in notes[:8])
    return (
        heading
        + '<p>Control pseudo-elements test the empirical null. At p &lt; 0.05, an approximately '
          '5% call rate is expected under a well-calibrated null. AUPRC/AUROC compare the same '
          'direct-target positives and non-targeting negatives for both methods.</p>'
        + '<h4>Empirical null</h4>' + data_table(calibration_rows, calibration_columns, limit=30)
        + '<h4>Direct-target discrimination</h4>' + data_table(auc_rows, auc_columns, limit=20)
        + (f'<ul>{note_text}</ul>' if note_text else "") + '</div>'
    )


def evaluation_artifact_content(root: Path | None) -> str:
    if not root:
        return ""
    evaluation = root / "evaluation_output"
    if not evaluation.is_dir():
        evaluation = root / "pipeline_dashboard" / "evaluation_output"
    if not evaluation.is_dir():
        return '<h3>Evaluation outputs</h3><div class="empty">Evaluation artifacts are not available yet.</div>'
    rows = []
    skip_notes = []
    for path in sorted(evaluation.iterdir()):
        if not path.is_file():
            continue
        if path.suffix == ".txt" and "skipped" in path.name:
            note = sanitized_tail(path, 20)
            if note:
                skip_notes.append(f'<div class="qc-callout"><strong>{html.escape(path.name)}</strong><pre>{html.escape(note)}</pre></div>')
        kind = path.suffix.lstrip(".").upper() or "file"
        rows.append({"artifact": path.name, "kind": kind, "bytes": path.stat().st_size})
    return (
        '<h3>Evaluation outputs and benchmark status</h3>' + "".join(skip_notes)
        + data_table(rows, [("artifact", "Artifact"), ("kind", "Type"), ("bytes", "Bytes")])
    )


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
    has_pipeline_error = bool(re.search(
        r"(?im)^(?:ERROR\s*~|.*Error executing process|.*Execution cancelled|"
        r".*terminated with an error|.*Session aborted\s*--\s*Cause)",
        log_tail,
    ))
    if not failed and not has_pipeline_error:
        return ""
    blocks = []
    for index, row in enumerate(rows):
        if row.get("status", "").upper() not in {"FAILED", "ABORTED"}:
            continue
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
        later_success = next((
            other for other in rows[index + 1:]
            if row.get("name") and other.get("name") == row["name"]
            and other.get("status", "").upper() in {"COMPLETED", "CACHED"}
        ), None)
        if later_success is not None:
            label = html.escape(later_success["status"].upper())
            runtime = html.escape(later_success.get("realtime") or later_success.get("duration") or "Unavailable")
            note += f'<p>Later trace entry with the same task name: <strong>{label}</strong> · runtime: {runtime}.</p>'
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
        memory = (row.get("peak_rss") or "").strip()
        if memory in {"", "-", "—"}:
            memory = "Unavailable"
        body.append(
            "<tr><td>" + html.escape(row.get("process", "").split(":")[-1]) + "</td>"
            "<td>" + html.escape(row.get("name", "")) + "</td>"
            f'<td><span class="pill {status}">{html.escape(row.get("status", "UNKNOWN"))}</span></td>'
            "<td>" + html.escape(row.get("realtime") or row.get("duration", "—")) + "</td>"
            "<td>" + html.escape(memory) + "</td></tr>"
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
        table = (
            '<h3>Observed guide capture by read pair</h3><p>One row is shown per sampled read pair. '
            'The row uses the first FASTQ name and the highest-scoring observed mate and orientation. '
            'Compare the observed configuration with the coordinates declared in SeqSpec.</p>' + table +
            '<details class="supporting-details"><summary>How to read the library structure QC</summary>'
            '<p>Total hits counts guide matches. Hit ratio divides that count by sampled reads for the candidate. '
            'Position purity is the fraction of hits at the most common start. Flank purity is the fraction of upstream '
            'sequences at that start with the most common flank. The score weights these terms and guide diversity.</p>'
            '<p>Look for a concentrated guide start, a consistent flank and guides spread across the reference. '
            'A mate carrying the cell barcode or UMI can have no guide hits. Few or dispersed hits on the expected '
            'guide mate merit checking the supplied library structure.</p></details>'
        )
    image = ""
    if image_path and image_path.exists() and image_path.stat().st_size <= 4_000_000:
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        image = f'<figure><img src="data:image/png;base64,{encoded}" alt="SeqSpec QC"><figcaption>SeqSpec read-structure QC</figcaption></figure>'
    return table + image


def continuation_content(report: dict[str, Any]) -> str:
    if not report:
        return ""
    parent = html.escape(str(report["parent_run_name"]))
    if report["qc_snapshot"]:
        provenance = "Mapping and QC below describe this checkpoint."
    elif report["qc_artifacts"]:
        provenance = "This older checkpoint has partial published QC. Some mapping metrics or plots may be unavailable."
    else:
        provenance = "Published mapping and QC files for this checkpoint are unavailable."
    cards = "".join(
        metric_card(label, display_value(report.get(key)), detail)
        for label, key, detail in (
            ("Checkpoint cells", "parent_cells", "after pipeline QC"),
            ("Submitted cells", "submitted_cells", "before the external allowlist"),
            ("Inference input cells", "accepted_cells", "after external selection"),
        )
    )
    covariates = ", ".join(html.escape(str(name)) for name in report.get("external_covariates", []))
    note = f"<p>External covariates: {covariates}</p>" if covariates else ""
    return (
        f'<section class="qc-callout"><h2>Continuation from {parent}</h2><p>{provenance}</p>'
        f'<div class="metrics">{cards}</div>{note}</section>'
    )


def render(args: argparse.Namespace) -> str:
    trace_rows = merge_live_tasks(
        read_trace(args.trace), getattr(args, "nextflow_log", None)
    )
    family_states = {family: family_state(trace_rows, family, args.status) for family, _, _ in FAMILIES}
    counts = Counter(row.get("status", "UNKNOWN").upper() for row in trace_rows)
    total_runtime = sum(duration_seconds(row.get("realtime") or row.get("duration", "")) for row in trace_rows)
    artifact_dir = getattr(args, "artifact_dir", None)
    checkpoint = artifact_dir / "checkpoint_dashboard" if artifact_dir else None
    continuation_record = checkpoint / "continuation.json" if checkpoint else None
    continuation = {}
    qc_metrics_json = getattr(args, "qc_metrics_json", None)
    guide_report = args.guide_report
    seqspec_table = args.seqspec_table
    if continuation_record and continuation_record.is_file():
        continuation = json.loads(continuation_record.read_text(encoding="utf-8"))
        qc_metrics_json = checkpoint / "pipeline_qc_metrics.json"
        guide_report = checkpoint / "guide_metadata.validation.json"
        seqspec_table = checkpoint / "guide_position_table.csv"
        for family in CHECKPOINT_FAMILIES:
            state = family_states[family]
            if not state["rows"] and continuation["qc_artifacts"]:
                state["status"] = "checkpoint"
    else:
        checkpoint = None
    checkpoint_finished = getattr(args, "checkpoint_ready", False) or any(
        row.get("process", "").split(":")[-1] == "preinference_checkpoint"
        and row.get("status", "").upper() in {"COMPLETED", "CACHED"}
        for row in trace_rows
    )
    ready = checkpoint_finished and args.status.lower() == "completed"
    final_preview = getattr(args, "final_dashboard_html", None)
    if continuation and args.status.lower() == "completed" and final_preview and final_preview.is_file():
        for family in ("inference", "evaluation", "final"):
            if family_states[family]["status"] == "pending":
                family_states[family]["status"] = "completed"
    if ready:
        for family in ("inference", "evaluation"):
            family_states[family]["status"] = "awaiting continuation"

    guide = json.loads(guide_report.read_text(encoding="utf-8")) if guide_report.is_file() else {}
    qc_data = {}
    if qc_metrics_json and qc_metrics_json.is_file():
        qc_data = json.loads(qc_metrics_json.read_text(encoding="utf-8"))
    elif continuation:
        observed = {}
        for key, filename in (("measurement_set_rna_qc", "measurement_set_qc_metrics.tsv"),
                              ("guide_mapping_qc", "guide_mapping_qc.tsv")):
            path = next(iter(sorted(checkpoint.rglob(filename))), None)
            if path:
                with path.open(newline="", encoding="utf-8") as table:
                    observed[key] = {"rows": list(csv.DictReader(table, delimiter="\t"))}
        observed["additional_qc"] = {"guide_mapping_qc": observed.pop("guide_mapping_qc", {})}
        qc_data = {
            "observed_metrics": observed,
            "dashboard_summary": {"filtering_summary": {"final_mudata": {"cells": continuation["parent_cells"]}}},
        }
    omitted_images = []
    artifact_images = collect_images(
        artifact_dir, getattr(args, "max_image_bytes", 32_000_000), omitted_images, checkpoint
    )
    guide_count = guide.get("row_count")
    guide_card = metric_card("Validated guide metadata rows", display_value(guide_count), "input validation report") if guide_count is not None else ""

    graph_nodes = []
    for index, (family, title, subtitle) in enumerate(FAMILIES, start=1):
        state = family_states[family]
        count_label = (
            state["status"] if (continuation and not state["rows"]) or state["status"] in {"checkpoint", "awaiting continuation"}
            else f'{state["completed"]} complete · {state["running"]} running · {state["failed"]} failed'
        )
        graph_nodes.append(
            f'<label for="family-tab-{family}" class="node {state["status"]}" data-family="{family}" tabindex="0" role="tab">'
            f'<span class="node-index">{index:02d}</span><span class="node-status"></span>'
            f'<strong>{html.escape(title)}</strong><small>{html.escape(subtitle)}</small>'
            f'<span class="node-count">{count_label}</span></label>'
        )

    sections = []
    for family, title, subtitle in FAMILIES:
        state = family_states[family]
        cards = [
            metric_card("Completed", state["completed"], "processes"),
            metric_card("Running", state["running"], "processes"),
            metric_card("Failed", state["failed"], "processes"),
            metric_card("Cached", state["cached"], "processes"),
            metric_card("Task runtime", fmt_seconds(state["runtime"]), "aggregate"),
        ]
        family_images = artifact_images.get(family, [])
        family_root = checkpoint if checkpoint and family in CHECKPOINT_FAMILIES else artifact_dir
        family_qc = {} if checkpoint and family in {"inference", "evaluation"} else qc_data
        extra = ""
        if family != "postconcat_qc":
            extra += qc_metrics_content(family_qc, family, family_root)
        if family == "mapping":
            mapping_cards = mapping_measurement_set_cards(qc_data)
            if mapping_cards:
                extra += '<details class="supporting-details"><summary>Mapping cards by measurement set</summary>' + mapping_cards + '</details>'
            extra += guide_mapping_qc_content(family_root)
            if any("sequencing_saturation" in path.stem for path, _ in family_images):
                extra += ('<h3>GEX sequencing saturation</h3><p>Saturation estimates 1 minus the fraction of usable reads '
                          'representing distinct barcode–UMI–gene molecules in RNA-QC-retained cells.</p>')
        if family == "preprocessing":
            measurement_cards, used_images = preprocessing_measurement_set_cards(qc_data, family_images)
            extra += measurement_cards
            family_images = [(path, encoded) for path, encoded in family_images if path not in used_images]
            extra += measurement_set_filter_flow_content(family_root)
        if family == "guide_assignment":
            assignment_cards, used_images = figure_measurement_set_cards(
                family_images, measurement_set_names(qc_data),
                "Guide assignment by measurement set",
                "Open a card to inspect guide-cell calling and assignment filtering for that set.",
                "guide calling and assignment filters",
            )
            extra += assignment_cards
            family_images = [(path, encoded) for path, encoded in family_images if path not in used_images]
            if not qc_data.get("observed_metrics", {}).get("additional_qc", {}).get("guide_assignment_filter", {}).get("rows"):
                extra += assignment_filter_flow_content(family_root, family)
        if family == 'postconcat_qc':
            process_images = [(path, encoded) for path, encoded in family_images if path.name == "postconcat_qc_flow.png"]
            extra += image_gallery(process_images, "Downstream filter flow")
            family_images = [(path, encoded) for path, encoded in family_images if path.name != "postconcat_qc_flow.png"]
            postconcat_cards, used_images = figure_measurement_set_cards(
                family_images, measurement_set_names(qc_data),
                "Post-concatenation views by measurement set",
                "Open a card to inspect the mitochondrial filter and available embedding and HTO plots for that set.",
                "measurement-set QC",
            )
            extra += postconcat_cards
            family_images = [(path, encoded) for path, encoded in family_images if path not in used_images]
            extra += assignment_filter_flow_content(family_root, family)
            extra += postconcat_qc_content(family_root)
            extra += qc_metrics_content(qc_data, family)
        if family == 'final':
            included = sum(len(items) for items in artifact_images.values())
            extra += f'<h3>QC image coverage</h3><p data-images-included="{included}" data-images-omitted="{len(omitted_images)}">{included} unique stage-specific QC images embedded; {len(omitted_images)} omitted by the configured size budget.</p>'
            if omitted_images:
                extra += data_table(omitted_images, [('plot', 'Omitted plot'), ('bytes', 'Bytes'), ('reason', 'Reason')], limit=len(omitted_images))
        if family == "input":
            if guide_card:
                extra = '<div class="metrics">' + guide_card + "".join(
                    metric_card(label, display_value(guide[key]), "guide metadata rows")
                    for key, label in (("targeting_rows", "Targeting"), ("control_rows", "Controls"))
                    if guide.get(key) is not None
                ) + '</div>' + extra
            else:
                extra = '<p>Guide metadata validation counts are unavailable for this report.</p>' + extra
        if family == "seqspec":
            seqspec_image = None if family_images or checkpoint else args.seqspec_image
            extra += seqspec_content(seqspec_table, seqspec_image)
        if family == "inference":
            extra += calibration_artifact_content(artifact_dir)
            extra += image_gallery(family_images, "Calibration and inference visualizations")
            family_images = []
            extra += ('<p>Inference is awaiting continuation from this checkpoint.</p>' if ready else
                      final_inference_content(final_preview))
        if family == "evaluation":
            extra += evaluation_artifact_content(artifact_dir)
        extra += image_gallery(family_images)
        task_details = (
            f'<div class="metrics telemetry">{"".join(cards)}</div><h3>Processes</h3>{process_table(state["rows"])}'
            if state["rows"] or (not continuation and not ready) else ""
        )
        sections.append(
            f'<section id="family-{family}" class="family-panel"><div class="family-heading">'
            f'<div><span class="eyebrow">Pipeline family</span><h2>{html.escape(title)}</h2>'
            f'<p>{html.escape(subtitle)}</p></div><span class="status-badge {state["status"]}">{state["status"]}</span></div>'
            f'{category_flow(family_qc, family)}{extra}{task_details}</section>'
        )

    if continuation and args.status.lower() == "completed":
        current = "inference"
    elif args.status.lower() == "completed" and measurement_set_names(qc_data):
        current = "preprocessing"
    else:
        current = next((family for family, _, _ in reversed(FAMILIES)
                        if family_states[family]["status"] in {"failed", "running", "completed"}),
                       "input")
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    family_tabs = "".join(
        f'<input class="family-tab" type="radio" name="pipeline-family" id="family-tab-{family}"'
        f'{" checked" if family == current else ""}>'
        for family, _, _ in FAMILIES
    )
    family_tab_css = "".join(
        f'#family-tab-{family}:checked~.graph-card .node[data-family="{family}"]'
        f'{{transform:translateY(-3px);border-color:var(--cyan);box-shadow:0 0 0 2px #0f8fc522}}'
        f'#family-tab-{family}:checked~.family-panels #family-{family}{{display:block}}'
        for family, _, _ in FAMILIES
    )
    failures = failure_content(
        trace_rows, getattr(args, "nextflow_log", None), getattr(args, "tail_lines", 30)
    )
    task_summary = "".join(
        metric_card(label, value, detail)
        for label, value, detail in (
            ("Completed", counts["COMPLETED"] + counts["CACHED"], "tasks"),
            ("Running", counts["RUNNING"] + counts["SUBMITTED"] + counts["NEW"], "tasks"),
            ("Failed", counts["FAILED"] + counts["ABORTED"], "tasks"),
            ("Cached", counts["CACHED"], "tasks"),
            ("Task runtime", fmt_seconds(total_runtime), "aggregate"),
        )
    ) if trace_rows or not continuation else ""
    document = f'''<!doctype html>
<html lang="en" data-dashboard="execution"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CRISPR Pipeline · {html.escape(args.run_name)}</title>
<style>
:root{{--bg:#f6f8fb;--panel:#ffffff;--panel2:#ffffff;--line:#d8e1ec;--text:#0f172a;--muted:#64748b;--cyan:#0f8fc5;--green:#16a36a;--red:#d1435b;--amber:#d68a00;--grey:#94a3b8}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--text);font:14px/1.5 Inter,ui-sans-serif,system-ui,sans-serif}}
.shell{{max-width:1500px;margin:auto;padding:28px}} header{{display:flex;justify-content:space-between;gap:24px;align-items:flex-start;margin-bottom:20px}} h1{{font-size:27px;margin:4px 0}} h2{{margin:3px 0 0;font-size:24px}} h3{{margin-top:28px}} p{{color:var(--muted);margin:4px 0}} .eyebrow{{color:var(--cyan);font:600 11px ui-monospace,monospace;letter-spacing:.14em;text-transform:uppercase}}
.run-id{{font-family:ui-monospace,monospace;color:var(--muted)}} .live{{display:flex;align-items:center;gap:8px;background:#fff;border:1px solid var(--line);border-radius:99px;padding:8px 12px}} .live i{{width:9px;height:9px;border-radius:50%;background:var(--amber);box-shadow:0 0 12px var(--amber)}} .live.completed i{{background:var(--green);box-shadow:0 0 12px var(--green)}} .live.failed i,.live.interrupted i{{background:var(--red);box-shadow:0 0 12px var(--red)}} .live.running i{{background:var(--cyan);box-shadow:0 0 12px var(--cyan)}}
header>div:first-child{{min-width:0}} header h1,.run-id{{overflow-wrap:anywhere}}
.execution-history summary,.failure-evidence summary{{overflow-wrap:anywhere}}
.summary,.metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:10px;margin:18px 0}} .metric{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px;box-shadow:0 5px 18px #0f172a0d}} .metric-label{{color:var(--muted);font-size:12px}} .metric-value{{font-size:24px;font-weight:750;margin-top:4px}} .metric-detail{{color:#66809d;font-size:11px}}
.telemetry{{margin:12px 0;gap:7px;grid-template-columns:repeat(auto-fit,minmax(115px,1fr))}} .telemetry .metric{{padding:7px 10px;box-shadow:none}} .telemetry .metric-value{{font-size:17px;margin-top:0}} .supporting-details{{border-color:var(--line);margin-bottom:14px}} .supporting-details>summary{{color:var(--text)}} a[data-open-family]{{color:#075985}}
.graph-card,.family-panel{{background:var(--panel);border:1px solid var(--line);border-radius:16px;padding:18px;margin-top:14px;box-shadow:0 12px 30px #0f172a12}} .graph-head{{display:flex;justify-content:space-between;align-items:center}} .graph{{display:flex;align-items:stretch;overflow-x:auto;padding:20px 2px 10px}} .node{{position:relative;text-decoration:none;flex:0 0 145px;min-height:132px;text-align:left;color:var(--text);background:#f8fafc;border:1px solid var(--line);border-radius:12px;padding:14px;cursor:pointer;transition:.18s}} .node:hover,.node.active{{transform:translateY(-3px);border-color:var(--cyan);box-shadow:0 0 0 2px #0f8fc522}} .node:not(:last-child){{margin-right:29px}} .node:not(:last-child):after{{content:'→';position:absolute;right:-23px;top:49px;color:#94a3b8;font-size:22px}} .node strong,.node small,.node-count{{display:block}} .node strong{{margin-top:18px}} .node small{{color:var(--muted);font-size:11px;min-height:34px}} .node-count{{font-size:10px;color:#7890aa;margin-top:7px}} .node-index{{font:600 10px ui-monospace,monospace;color:#6c86a1}} .node-status{{position:absolute;right:12px;top:12px;width:10px;height:10px;border-radius:50%;background:var(--grey)}}
.node.completed .node-status,.completed.status-badge,.node.checkpoint .node-status,.checkpoint.status-badge{{background:var(--green)}} .node.running .node-status,.running.status-badge{{background:var(--cyan);box-shadow:0 0 12px var(--cyan)}} .node.failed .node-status,.failed.status-badge{{background:var(--red)}} .node.pending{{opacity:.65}} .family-tab{{position:absolute;inline-size:1px;block-size:1px;opacity:0;pointer-events:none}} .family-panel{{display:none;scroll-margin-top:12px}} {family_tab_css} .family-heading{{display:flex;justify-content:space-between;align-items:flex-start}} .status-badge{{border-radius:99px;padding:5px 10px;text-transform:uppercase;font-size:10px;font-weight:800;color:#06121e;background:var(--grey)}}
.table-wrap{{overflow:auto;border:1px solid var(--line);border-radius:10px}} table{{border-collapse:collapse;width:100%;min-width:700px}} th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid #e2e8f0}} th{{color:#475569;background:#f1f5f9;font-size:11px;text-transform:uppercase;letter-spacing:.06em}} td{{font-family:ui-monospace,monospace;font-size:12px}} .pill{{padding:3px 7px;border-radius:99px;background:#e2e8f0;font-size:10px}} .pill.completed,.pill.cached{{background:#dcfce7;color:#166534}} .pill.running,.pill.submitted,.pill.new{{background:#e0f2fe;color:#075985}} .pill.failed,.pill.aborted{{background:#ffe4e6;color:#be123c}} figure{{margin:18px 0;background:#fff;border-radius:12px;padding:10px}} figure img{{display:block;max-width:100%;margin:auto}} figcaption{{color:#50647b;padding:8px 4px 2px}} .empty{{color:var(--muted);border:1px dashed var(--line);border-radius:10px;padding:20px}} footer{{color:#607994;font-size:11px;margin:20px 2px}}
.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}} .gallery figure{{margin:0;min-width:0}} .failure-evidence{{background:#fff1f2;border:1px solid #fecdd3;border-radius:16px;padding:18px;margin-top:14px;box-shadow:0 12px 30px #0f172a12}} details{{background:#fff;border:1px solid #fecdd3;border-radius:10px;margin-top:10px;padding:10px 12px}} summary{{cursor:pointer;font-weight:700;color:#be123c}} pre{{white-space:pre-wrap;word-break:break-word;max-height:340px;overflow:auto;background:#f8fafc;border-radius:8px;padding:12px;color:#334155;font:11px/1.45 ui-monospace,monospace}} .evidence-note{{color:var(--amber)}} .table-search{{width:min(520px,100%);background:#fff;color:var(--text);border:1px solid var(--line);border-radius:9px;padding:10px 12px;margin:0 0 10px}} .result-block{{border-top:1px solid var(--line);margin-top:22px;padding-top:2px}} .qc-callout{{border:1px solid #fde68a;background:#fffbeb;border-radius:10px;padding:12px;margin:10px 0}}
.measurement-section{{margin:24px 0;padding:16px;border:1px solid var(--line);border-radius:14px;background:#f8fafc}} .measurement-section-head{{display:flex;align-items:flex-start;justify-content:space-between;gap:18px;margin-bottom:10px}} .measurement-section-head h3{{margin:3px 0}} .measurement-count{{white-space:nowrap;background:#dff3fb;color:#075985;border-radius:99px;padding:6px 10px;font-size:11px;font-weight:800}}
details.measurement-card{{border-color:var(--line);border-radius:12px;padding:0;overflow:hidden;box-shadow:0 4px 14px #0f172a0a}} details.measurement-card[open]{{border-color:#7dd3fc;box-shadow:0 0 0 2px #0ea5e91a}} details.measurement-card>summary{{list-style:none;color:var(--text);display:grid;grid-template-columns:auto minmax(240px,1fr) auto auto;align-items:center;gap:18px;padding:14px 16px;background:#fff}} details.measurement-card>summary::-webkit-details-marker{{display:none}} details.measurement-card>summary:before{{content:'›';font-size:22px;color:var(--cyan);transition:.15s}} details.measurement-card[open]>summary:before{{transform:rotate(90deg)}}
.measurement-name,.measurement-subtitle,.measurement-summary span{{display:block}} .measurement-name{{font:700 13px ui-monospace,monospace}} .measurement-subtitle{{color:var(--muted);font-size:11px;font-weight:500;margin-top:2px}} .measurement-summary{{text-align:right}} .measurement-summary strong{{font-size:18px}} .measurement-summary span{{color:var(--muted);font-size:10px}} .modality-chips{{display:flex;gap:5px;flex-wrap:wrap;justify-content:flex-end}} .modality-chips span{{background:#e0f2fe;color:#075985;border-radius:99px;padding:4px 7px;font-size:10px;text-transform:uppercase}}
.measurement-body{{border-top:1px solid var(--line);padding:16px;background:#fbfdff}} .measurement-body h4{{margin:4px 0 10px;text-transform:capitalize}} .modality-block{{border:1px solid var(--line);border-radius:10px;background:#fff;padding:14px;margin-bottom:12px}} .measurement-facts{{display:grid;grid-template-columns:repeat(auto-fit,minmax(125px,1fr));gap:8px;margin:10px 0}} .measurement-fact{{background:#f1f5f9;border-radius:8px;padding:9px}} .measurement-fact span,.measurement-fact strong{{display:block}} .measurement-fact span{{color:var(--muted);font-size:10px}} .measurement-fact strong{{font-size:13px;margin-top:3px}}
.mini-bar{{margin:10px 0}} .mini-bar-label{{display:flex;justify-content:space-between;font-size:11px;color:var(--muted)}} .mini-bar-track{{height:7px;border-radius:99px;background:#e2e8f0;overflow:hidden;margin-top:4px}} .mini-bar-track i{{display:block;height:100%;background:linear-gradient(90deg,#22d3ee,#0f8fc5);border-radius:99px}} .filter-stages{{display:flex;gap:22px;overflow-x:auto;padding:3px 2px 12px}} .filter-stage{{position:relative;flex:1 0 145px;border:1px solid var(--line);border-radius:9px;background:#fff;padding:10px}} .filter-stage:not(:last-child):after{{content:'→';position:absolute;right:-17px;top:31%;color:#94a3b8}} .filter-stage span,.filter-stage strong{{display:block}} .filter-stage span{{font-size:10px;color:var(--muted)}} .filter-stage strong{{font-size:16px;margin-top:3px}}
.measurement-plots{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:12px;margin-top:14px}} .measurement-plots figure{{margin:0;border:1px solid var(--line)}}
.process-flow{{border:1px solid var(--line);background:#f8fafc;border-radius:12px;padding:14px;margin:18px 0}} .flow-steps{{display:flex;align-items:stretch;overflow-x:auto;gap:24px;padding:4px 2px}} .flow-step{{position:relative;flex:1 0 180px;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px}} .flow-step:not(:last-child):after{{content:'→';position:absolute;right:-18px;top:38%;color:#94a3b8;font-size:20px}} .flow-step strong,.flow-step span{{display:block}} .flow-step span{{color:var(--muted);font-size:11px;margin-top:5px}} .filter-chips{{display:flex;flex-wrap:wrap;gap:7px}} .filter-chips span{{background:#e0f2fe;color:#0c4a6e;border-radius:99px;padding:5px 9px;font-size:11px}} .flow-note{{margin-top:9px;font-size:12px}}
@media(max-width:700px){{.shell{{padding:15px}}header{{display:block}}.live{{margin-top:12px;width:max-content}}details.measurement-card>summary{{grid-template-columns:auto 1fr}}.measurement-summary,.modality-chips{{grid-column:2;text-align:left;justify-content:flex-start}}.measurement-plots{{grid-template-columns:1fr}}}}
</style></head><body><div class="shell">
<header><div><span class="eyebrow">CRISPR Pipeline · execution dashboard</span><h1>{html.escape(args.run_name)}</h1><div class="run-id">{html.escape(args.run_id)}</div></div><div class="live {html.escape(args.status.lower())}"><i></i><span>{"READY FOR INFERENCE" if ready else html.escape(args.status.upper())}</span></div></header>
{continuation_content(continuation)}
{family_tabs}
{f'<div class="summary telemetry">{task_summary}{guide_card}</div>' if task_summary or guide_card else ""}
<div class="graph-card"><div class="graph-head"><div><span class="eyebrow">Live dependency view</span><h2>Pipeline execution steps</h2></div><p>Select a stage to inspect its QC and tasks</p></div><nav class="graph">{"".join(graph_nodes)}</nav></div>
{failures}<div class="family-panels">{"".join(sections)}</div><footer>Generated {generated} · Self-contained W&amp;B HTML media · No credentials, FASTQs or unbounded task logs embedded</footer></div>
<script>
document.querySelectorAll("[data-open-family]").forEach(link => {{
    link.addEventListener("click", () => {{ document.getElementById("family-tab-" + link.dataset.openFamily).checked = true; }});
}});
document.querySelectorAll(".node").forEach(node => {{
    node.addEventListener("keydown", event => {{
        if (event.key === "Enter" || event.key === " ") {{ event.preventDefault(); node.click(); }}
    }});
}});
document.querySelectorAll(".table-search").forEach(input => {{
    const table = document.getElementById(input.getAttribute("aria-controls"));
    input.addEventListener("input", () => {{
        const query = input.value.trim().toLowerCase();
        for (const row of table.tBodies[0].rows) row.hidden = !row.textContent.toLowerCase().includes(query);
    }});
}});
</script>
</body></html>'''
    return enhance_dashboard_html(document, args.status)


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
    parser.add_argument("--final-dashboard-html", type=Path)
    parser.add_argument("--nextflow-log", type=Path)
    parser.add_argument("--tail-lines", type=int, default=30)
    parser.add_argument("--max-image-bytes", type=int, default=38_000_000)
    parser.add_argument("--max-html-bytes", type=int, default=100_000_000)
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
