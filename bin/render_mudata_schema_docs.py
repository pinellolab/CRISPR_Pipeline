#!/usr/bin/env python3
"""Render docs/mudata_schema.md from the canonical MuData field catalog."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "docs" / "mudata_schema_catalog.tsv"
OUTPUT = ROOT / "docs" / "mudata_schema.md"
GOOGLE_SHEET_URL = (
    "https://docs.google.com/spreadsheets/d/"
    "1hwGyxCtwwgnpzEgdt7BlkK1tuJ22BPMpY-BEtCpTKm0/edit"
)


def esc(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def main() -> None:
    with CATALOG.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sections: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        sections.setdefault(row["sheet"], []).append(row)

    lines = [
        "# MuData field reference",
        "",
        "The final `pipeline_outputs/inference_mudata.h5mu` is a MuData object whose rows are retained cell barcodes and whose principal modalities are `gene` and `guide`. A third `hashing` modality is present when data hashing is enabled. RNA and guide `.X` matrices contain raw UMI counts; the binary guide calls used for inference are stored separately in `guide.layers['guide_assignment']`.",
        "",
        f"The canonical multi-tab field dictionary is also available as a [Google Sheet]({GOOGLE_SHEET_URL}). The checked-in source for both representations is [`mudata_schema_catalog.tsv`](mudata_schema_catalog.tsv). Regenerate this page with `python3 bin/render_mudata_schema_docs.py`.",
        "",
        "## Scope and stability",
        "",
        "Fields marked **Always** are part of the current pipeline contract. Optional fields depend on hashing, Scrublet, guide-assignment method, inference method, or available metadata. The pipeline deliberately preserves additional `guide.var` columns from the input guide metadata; therefore dataset-specific columns may appear beyond the currently observed extensions documented below. Current output keys use `local_analysis_*` and `global_analysis_*`; `cis_*` and `trans_*` are legacy aliases retained for interpretation of older files.",
        "",
        "The catalog was audited against `dev` source commit `33d068dad7f8163c313bd6fa62aa45c6b35bffa2` and completed SCEPTRE, CLEANSER, enhancer-screen, Gasperini, and Replogle outputs.",
        "",
    ]
    for section, group in sections.items():
        lines += [f"## {section}", "", "| Path | Field | Type | Availability | Description | Producer/source |", "|---|---|---|---|---|---|"]
        for row in group:
            lines.append("| " + " | ".join(esc(row[key]) for key in ("path", "field", "type", "availability", "description", "producer_or_source")) + " |")
        lines.append("")
    OUTPUT.write_text("\n".join(lines))
    print(f"Wrote {OUTPUT} with {len(rows)} catalog rows")


if __name__ == "__main__":
    main()
