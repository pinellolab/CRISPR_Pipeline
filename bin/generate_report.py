#!/usr/bin/env python3
"""
Generate an HTML/PDF summary report of QC results.

This script combines metrics tables and plots from the QC pipeline into
a single HTML document that can optionally be converted to PDF.

Inputs expected:
  - mapping_gene_dir: gene_metrics.tsv, gene_histograms_by_batch.png, etc.
  - mapping_guide_dir: guide_metrics.tsv, guide_histograms_by_batch.png, etc.
  - intended_target_dir: intended_target_metrics.tsv, intended_target_volcano.png, etc.

Output:
  - HTML report (and optionally PDF via WeasyPrint)
"""

import argparse
import base64
import logging
import os
import sys
from datetime import datetime
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


# ---------------------------------------------------------------------
# CSS Styling
# ---------------------------------------------------------------------
CSS = """
<style>
    @page {
        size: letter portrait;
        margin: 0.75in;
    }

    body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        font-size: 11pt;
        line-height: 1.5;
        color: #333;
        max-width: 8in;
        margin: 0 auto;
        padding: 20px;
    }

    h1 {
        font-size: 28pt;
        font-weight: 600;
        color: #2c3e50;
        text-align: center;
        margin-bottom: 5px;
        page-break-after: avoid;
    }

    h2 {
        font-size: 18pt;
        font-weight: 600;
        color: #2c3e50;
        border-bottom: 2px solid #3498db;
        padding-bottom: 8px;
        margin-top: 40px;
        margin-bottom: 15px;
        page-break-after: avoid;
    }

    h3 {
        font-size: 14pt;
        font-weight: 600;
        color: #34495e;
        margin-top: 25px;
        margin-bottom: 10px;
        page-break-after: avoid;
    }

    .subtitle {
        font-size: 14pt;
        color: #666;
        text-align: center;
        margin-bottom: 5px;
    }

    .date {
        font-size: 11pt;
        color: #888;
        text-align: center;
        margin-bottom: 30px;
    }

    .description {
        font-size: 11pt;
        color: #555;
        font-style: italic;
        margin-bottom: 15px;
        background-color: #f8f9fa;
        padding: 10px 15px;
        border-radius: 5px;
        border-left: 3px solid #3498db;
    }

    table {
        width: 100%;
        border-collapse: collapse;
        margin: 15px 0 25px 0;
        font-size: 10pt;
        page-break-inside: avoid;
    }

    table th {
        background-color: #2c3e50;
        color: white;
        font-weight: 600;
        padding: 10px 8px;
        text-align: center;
        border: 1px solid #2c3e50;
    }

    table td {
        padding: 8px;
        text-align: center;
        border: 1px solid #ddd;
    }

    table tr:nth-child(even) {
        background-color: #f8f9fa;
    }

    table tr:hover {
        background-color: #e8f4f8;
    }

    .figure-container {
        text-align: center;
        margin: 20px 0;
        page-break-inside: avoid;
    }

    .figure-container img {
        max-width: 100%;
        height: auto;
        border: 1px solid #ddd;
        border-radius: 4px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }

    .figure-caption {
        font-size: 10pt;
        color: #666;
        font-style: italic;
        margin-top: 8px;
    }

    .two-column {
        display: flex;
        justify-content: space-between;
        gap: 20px;
        margin: 20px 0;
        page-break-inside: avoid;
    }

    .two-column .figure-container {
        flex: 1;
        margin: 0;
    }

    .two-column img {
        max-width: 100%;
    }

    hr {
        border: none;
        border-top: 1px solid #ddd;
        margin: 30px 0;
    }

    .page-break {
        page-break-before: always;
    }
</style>
"""


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def image_to_base64(image_path: str) -> Optional[str]:
    """Convert an image file to base64 string for HTML embedding."""
    if not os.path.exists(image_path):
        logger.warning(f"Image not found: {image_path}")
        return None

    with open(image_path, "rb") as f:
        data = base64.b64encode(f.read()).decode("utf-8")

    ext = os.path.splitext(image_path)[1].lower()
    mime_types = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".gif": "image/gif"}
    mime = mime_types.get(ext, "image/png")

    return f"data:{mime};base64,{data}"


def format_dataframe(
    df: pd.DataFrame,
    column_renames: Optional[dict] = None,
    column_formats: Optional[dict] = None,
) -> pd.DataFrame:
    """Format a DataFrame for HTML display."""
    df = df.copy()

    # Apply number formatting
    if column_formats:
        for col, fmt in column_formats.items():
            if col in df.columns:
                if fmt == "pct":
                    df[col] = df[col].apply(lambda x: f"{x*100:.1f}%" if pd.notna(x) else "-")
                elif fmt == "int":
                    df[col] = df[col].apply(lambda x: f"{int(x):,}" if pd.notna(x) else "-")
                elif fmt == "float1":
                    df[col] = df[col].apply(lambda x: f"{x:.1f}" if pd.notna(x) else "-")
                elif fmt == "float2":
                    df[col] = df[col].apply(lambda x: f"{x:.2f}" if pd.notna(x) else "-")
                elif fmt == "float3":
                    df[col] = df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "-")

    # Rename columns
    if column_renames:
        df = df.rename(columns=column_renames)

    return df


def make_table_html(
    tsv_path: str,
    columns_to_show: Optional[list] = None,
    column_renames: Optional[dict] = None,
    column_formats: Optional[dict] = None,
) -> str:
    """Load a TSV and convert to HTML table."""
    if not os.path.exists(tsv_path):
        logger.warning(f"File not found: {tsv_path}")
        return "<p><em>Data not available</em></p>"

    df = pd.read_csv(tsv_path, sep="\t")

    # Select columns
    if columns_to_show:
        df = df[[c for c in columns_to_show if c in df.columns]]

    # Format
    df = format_dataframe(df, column_renames=column_renames, column_formats=column_formats)

    return df.to_html(index=False, border=0, classes="metrics-table")


def make_figure_html(image_path: str, caption: str = "") -> str:
    """Create HTML for a figure with optional caption."""
    b64 = image_to_base64(image_path)
    if b64 is None:
        return "<p><em>Figure not available</em></p>"

    html = f'<div class="figure-container">\n'
    html += f'  <img src="{b64}" alt="{caption}">\n'
    if caption:
        html += f'  <div class="figure-caption">{caption}</div>\n'
    html += '</div>\n'
    return html


def make_two_figures_html(
    image_path1: str,
    image_path2: str,
    caption1: str = "",
    caption2: str = "",
) -> str:
    """Create HTML for two figures side by side."""
    b64_1 = image_to_base64(image_path1)
    b64_2 = image_to_base64(image_path2)

    html = '<div class="two-column">\n'

    # Left figure
    html += '  <div class="figure-container">\n'
    if b64_1:
        html += f'    <img src="{b64_1}" alt="{caption1}">\n'
        if caption1:
            html += f'    <div class="figure-caption">{caption1}</div>\n'
    else:
        html += '    <p><em>Figure not available</em></p>\n'
    html += '  </div>\n'

    # Right figure
    html += '  <div class="figure-container">\n'
    if b64_2:
        html += f'    <img src="{b64_2}" alt="{caption2}">\n'
        if caption2:
            html += f'    <div class="figure-caption">{caption2}</div>\n'
    else:
        html += '    <p><em>Figure not available</em></p>\n'
    html += '  </div>\n'

    html += '</div>\n'
    return html


# ---------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------
def generate_qc_report(
    output_path: str,
    sample_name: str = "Sample",
    mapping_gene_dir: Optional[str] = None,
    mapping_guide_dir: Optional[str] = None,
    intended_target_dir: Optional[str] = None,
    generate_pdf: bool = False,
) -> None:
    """
    Generate an HTML QC summary report, optionally converting to PDF.

    Parameters
    ----------
    output_path : str
        Path for output file (.html or .pdf).
    sample_name : str
        Sample name for report title.
    mapping_gene_dir : str, optional
        Directory containing gene mapping QC results.
    mapping_guide_dir : str, optional
        Directory containing guide mapping QC results.
    intended_target_dir : str, optional
        Directory containing intended target QC results.
    generate_pdf : bool
        If True, also generate PDF using WeasyPrint.
    """
    logger.info(f"Generating QC report: {output_path}")

    timestamp = datetime.now().strftime("%B %d, %Y")

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QC Report - {sample_name}</title>
    {CSS}
</head>
<body>

<h1>QC Summary Report</h1>
<div class="subtitle">{sample_name}</div>
<div class="date">Generated: {timestamp}</div>

<hr>
"""

    # -----------------------------------------------------------------
    # Section 1: Gene Expression Mapping QC
    # -----------------------------------------------------------------
    if mapping_gene_dir and os.path.isdir(mapping_gene_dir):
        html += """
<h2>1. Gene Expression Mapping QC</h2>

<div class="description">
Quality metrics for gene expression capture including UMI counts per cell,
number of detected genes, and mitochondrial content. These metrics help identify
low-quality cells and batch effects.
</div>

<h3>Summary Statistics</h3>
"""
        html += make_table_html(
            os.path.join(mapping_gene_dir, "gene_metrics.tsv"),
            columns_to_show=["batch", "n_cells", "umi_median", "umi_mean", "genes_median", "genes_mean", "mito_median", "mito_mean"],
            column_renames={
                "batch": "Batch",
                "n_cells": "Cells",
                "umi_median": "UMI (Median)",
                "umi_mean": "UMI (Mean)",
                "genes_median": "Genes (Median)",
                "genes_mean": "Genes (Mean)",
                "mito_median": "Mito % (Median)",
                "mito_mean": "Mito % (Mean)",
            },
            column_formats={
                "n_cells": "int",
                "umi_median": "int",
                "umi_mean": "int",
                "genes_median": "int",
                "genes_mean": "int",
                "mito_median": "float2",
                "mito_mean": "float2",
            },
        )

        html += "\n<h3>Distributions by Batch</h3>\n"
        html += make_figure_html(
            os.path.join(mapping_gene_dir, "gene_histograms_by_batch.png"),
            "Distribution of UMI counts, detected genes, and mitochondrial percentage per cell, stratified by batch."
        )

        html += "\n<h3>Cell Counts and UMI Distribution</h3>\n"
        html += make_two_figures_html(
            os.path.join(mapping_gene_dir, "gene_cells_per_batch.png"),
            os.path.join(mapping_gene_dir, "gene_knee_plot.png"),
            "Number of cells per batch",
            "Knee plot of UMI counts (ranked barcodes)"
        )

    # -----------------------------------------------------------------
    # Section 2: Guide Mapping QC
    # -----------------------------------------------------------------
    if mapping_guide_dir and os.path.isdir(mapping_guide_dir):
        html += """
<div class="page-break"></div>
<h2>2. Guide Mapping QC</h2>

<div class="description">
Quality metrics for guide RNA capture and assignment. This section evaluates guide UMI counts
per cell and the number of guides assigned to each cell, which affects the interpretability
of perturbation effects.
</div>

<h3>Summary Statistics</h3>
"""
        html += make_table_html(
            os.path.join(mapping_guide_dir, "guide_metrics.tsv"),
            columns_to_show=["batch", "n_cells", "guide_umi_median", "n_cells_with_guide", "n_cells_exactly_1_guide", "guides_per_cell_mean", "guides_per_cell_max"],
            column_renames={
                "batch": "Batch",
                "n_cells": "Cells",
                "guide_umi_median": "Guide UMI (Median)",
                "n_cells_with_guide": "Cells w/ Guide",
                "n_cells_exactly_1_guide": "Cells w/ 1 Guide",
                "guides_per_cell_mean": "Guides/Cell (Mean)",
                "guides_per_cell_max": "Guides/Cell (Max)",
            },
            column_formats={
                "n_cells": "int",
                "guide_umi_median": "int",
                "n_cells_with_guide": "int",
                "n_cells_exactly_1_guide": "int",
                "guides_per_cell_mean": "float2",
            },
        )

        html += "\n<h3>Distributions by Batch</h3>\n"
        html += make_figure_html(
            os.path.join(mapping_guide_dir, "guide_histograms_by_batch.png"),
            "Distribution of guide UMI counts and number of guides assigned per cell, stratified by batch."
        )

        html += "\n<h3>Guide Assignment Overview</h3>\n"
        html += make_figure_html(
            os.path.join(mapping_guide_dir, "guide_histograms.png"),
            "Overall distributions including guide UMIs per cell, guides per cell, and cells per guide (colored by guide type)."
        )

    # -----------------------------------------------------------------
    # Section 3: Intended Target Inference QC
    # -----------------------------------------------------------------
    if intended_target_dir and os.path.isdir(intended_target_dir):
        html += """
<div class="page-break"></div>
<h2>3. Intended Target Inference QC</h2>

<div class="description">
Evaluation of guide knockdown efficiency on their intended target genes. This analysis uses
cis-regulatory differential expression testing to measure how effectively each guide reduces
expression of its target gene compared to non-targeting controls.
</div>

<h3>Knockdown Efficiency Summary</h3>
"""
        html += make_table_html(
            os.path.join(intended_target_dir, "intended_target_metrics.tsv"),
            columns_to_show=["n_guides_tested", "n_strong_knockdowns", "frac_strong_knockdowns", "n_significant", "frac_significant", "median_log2fc", "mean_log2fc"],
            column_renames={
                "n_guides_tested": "Guides Tested",
                "n_strong_knockdowns": "Strong KD (≥60%)",
                "frac_strong_knockdowns": "% Strong KD",
                "n_significant": "Significant (p<0.05)",
                "frac_significant": "% Significant",
                "median_log2fc": "Median Log2FC",
                "mean_log2fc": "Mean Log2FC",
            },
            column_formats={
                "n_guides_tested": "int",
                "n_strong_knockdowns": "int",
                "frac_strong_knockdowns": "pct",
                "n_significant": "int",
                "frac_significant": "pct",
                "median_log2fc": "float3",
                "mean_log2fc": "float3",
            },
        )

        html += "\n<h3>Knockdown Results</h3>\n"
        html += make_two_figures_html(
            os.path.join(intended_target_dir, "intended_target_volcano.png"),
            os.path.join(intended_target_dir, "intended_target_log2fc_distribution.png"),
            "Volcano plot: effect size (log2FC) vs. statistical significance (-log10 p-value)",
            "Distribution of log2 fold changes across all tested guides"
        )

    # Close HTML
    html += """
</body>
</html>
"""

    # Write HTML
    html_path = output_path if output_path.endswith(".html") else output_path.replace(".pdf", ".html")
    with open(html_path, "w") as f:
        f.write(html)
    logger.info(f"Saved HTML report to: {html_path}")

    # Generate PDF if requested
    if generate_pdf or output_path.endswith(".pdf"):
        try:
            from weasyprint import HTML
            pdf_path = output_path if output_path.endswith(".pdf") else output_path.replace(".html", ".pdf")
            HTML(string=html).write_pdf(pdf_path)
            logger.info(f"Saved PDF report to: {pdf_path}")
        except ImportError:
            logger.warning("WeasyPrint not installed. Run 'pip install weasyprint' to generate PDFs.")
        except Exception as e:
            logger.error(f"Failed to generate PDF: {e}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an HTML/PDF summary report of QC results."
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Path for output file (.html or .pdf)."
    )
    parser.add_argument(
        "--sample-name", default="Sample",
        help="Sample name for report title (default: 'Sample')."
    )
    parser.add_argument(
        "--mapping-gene-dir",
        help="Directory containing gene mapping QC results."
    )
    parser.add_argument(
        "--mapping-guide-dir",
        help="Directory containing guide mapping QC results."
    )
    parser.add_argument(
        "--intended-target-dir",
        help="Directory containing intended target QC results."
    )
    parser.add_argument(
        "--pdf", action="store_true",
        help="Also generate PDF output (requires WeasyPrint)."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    generate_qc_report(
        output_path=args.output,
        sample_name=args.sample_name,
        mapping_gene_dir=args.mapping_gene_dir,
        mapping_guide_dir=args.mapping_guide_dir,
        intended_target_dir=args.intended_target_dir,
        generate_pdf=args.pdf,
    )
