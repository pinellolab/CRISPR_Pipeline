#!/usr/bin/env python3
"""Generate 10x-style sequencing-saturation curves from corrected BUS files.

For each downsampling fraction p, a molecule represented by r usable reads is
expected to remain with probability 1-(1-p)^r.  Sequencing saturation is then
1 - expected unique (cell barcode, UMI, gene) molecules / expected usable
reads.  This is the same duplicate-read definition documented by 10x, applied
analytically rather than by materializing downsampled FASTQs.
"""

import argparse
import json
import re
import struct
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BUS_DTYPE = np.dtype(
    [
        ("barcode", "<u8"),
        ("umi", "<u8"),
        ("ec", "<i4"),
        ("count", "<u4"),
        ("flags", "<u4"),
        ("pad", "<u4"),
    ]
)


def encode_dna(sequence):
    value = 0
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    for base in str(sequence).upper():
        if base not in mapping:
            raise ValueError(f"Non-ACGT barcode cannot be encoded: {sequence}")
        value = (value << 2) | mapping[base]
    return value


def read_bus_header(handle):
    if handle.read(4) != b"BUS\x00":
        raise ValueError("Not a BUS file")
    version, barcode_length, umi_length, text_length = struct.unpack("<IIII", handle.read(16))
    text = handle.read(text_length).decode("utf-8", errors="replace")
    return version, barcode_length, umi_length, text


def transcript_to_gene(t2g_path):
    table = pd.read_csv(t2g_path, sep="\t", header=None, comment="#", dtype=str)
    if table.shape[1] < 2:
        raise ValueError("Transcript-to-gene file must have at least two columns")
    return dict(zip(table.iloc[:, 0], table.iloc[:, 1]))


def ec_gene_map(mapping_dir, t2g_path):
    tx_to_gene = transcript_to_gene(t2g_path)
    transcripts = (mapping_dir / "transcripts.txt").read_text().splitlines()
    transcript_genes = np.array([tx_to_gene.get(tx, "") for tx in transcripts], dtype=object)
    lines = (mapping_dir / "matrix.ec").read_text().splitlines()
    parsed = [(int(line.split("\t", 1)[0]), line.split("\t", 1)[1]) for line in lines]
    max_ec = max((ec_id for ec_id, _tx_text in parsed), default=-1)
    result = np.full(max_ec + 1, -1, dtype=np.int32)
    gene_ids = {}
    for ec_id, tx_text in parsed:
        genes = {transcript_genes[int(i)] for i in tx_text.split(",") if transcript_genes[int(i)]}
        if len(genes) == 1:
            gene = next(iter(genes))
            result[ec_id] = gene_ids.setdefault(gene, len(gene_ids))
    return result, len(gene_ids)


def batch_name(mapping_dir):
    match = re.match(r"(.+)_ks_transcripts_out$", mapping_dir.name)
    return match.group(1) if match else mapping_dir.name


def selected_barcodes(filtered, batch, covariates):
    if "batch" not in filtered.obs:
        if len(covariates) != 1:
            raise ValueError("Filtered AnnData lacks batch metadata for a multi-batch run")
        subset = filtered.obs
    else:
        subset = filtered.obs.loc[filtered.obs["batch"].astype(str) == str(batch)]
    if subset.empty:
        raise ValueError(f"No filtered cells found for mapping batch {batch}")
    if "corrected_barcode" in subset:
        raw = subset["corrected_barcode"].astype(str)
    else:
        rows = covariates.loc[covariates["batch"].astype(str) == str(batch)]
        if len(rows) != 1:
            raise ValueError(f"Expected one barcode_key for batch {batch}, found {len(rows)}")
        key = str(rows.iloc[0]["barcode_key"])
        suffix = "_" + key
        raw = pd.Index(subset.index.astype(str)).map(
            lambda name: name[: -len(suffix)] if name.endswith(suffix) else name
        )
    return {encode_dna(barcode) for barcode in raw}


def _survival(read_counts, fractions):
    values = np.asarray(read_counts, dtype=np.float64)
    if values.size == 0:
        return np.zeros(len(fractions), dtype=np.float64)
    return np.sum(1.0 - np.power(1.0 - fractions[:, None], values[None, :]), axis=1)


def analyze_bus(bus_path, selected, ec_to_gene, fractions, chunk_records=1_000_000):
    selected_sorted = np.array(sorted(selected), dtype=np.uint64)
    per_cell_umis = []
    per_cell_genes = []
    total_reads = 0.0
    total_expected_umis = np.zeros(len(fractions), dtype=np.float64)
    current_barcode = None
    current_umi = None
    umi_gene_reads = {}
    cell_molecule_reads = []
    cell_gene_reads = {}

    def flush_umi():
        nonlocal umi_gene_reads, cell_molecule_reads, cell_gene_reads
        for gene, count in umi_gene_reads.items():
            cell_molecule_reads.append(count)
            cell_gene_reads[gene] = cell_gene_reads.get(gene, 0) + count
        umi_gene_reads = {}

    def flush_cell():
        nonlocal cell_molecule_reads, cell_gene_reads, total_reads, total_expected_umis
        if not cell_molecule_reads:
            return
        molecule_curve = _survival(cell_molecule_reads, fractions)
        gene_curve = _survival(list(cell_gene_reads.values()), fractions)
        per_cell_umis.append(molecule_curve)
        per_cell_genes.append(gene_curve)
        total_expected_umis += molecule_curve
        total_reads += float(np.sum(cell_molecule_reads))
        cell_molecule_reads = []
        cell_gene_reads = {}

    with open(bus_path, "rb") as handle:
        _version, _barcode_length, _umi_length, _text = read_bus_header(handle)
        while True:
            records = np.fromfile(handle, dtype=BUS_DTYPE, count=chunk_records)
            if records.size == 0:
                break
            positions = np.searchsorted(selected_sorted, records["barcode"])
            valid_position = positions < len(selected_sorted)
            selected_mask = np.zeros(records.size, dtype=bool)
            selected_mask[valid_position] = (
                selected_sorted[positions[valid_position]] == records["barcode"][valid_position]
            )
            records = records[selected_mask]
            if records.size == 0:
                continue
            good_ec = (records["ec"] >= 0) & (records["ec"] < len(ec_to_gene))
            records = records[good_ec]
            if records.size == 0:
                continue
            genes = ec_to_gene[records["ec"]]
            usable = genes >= 0
            records = records[usable]
            genes = genes[usable]
            for barcode, umi, gene, count in zip(
                records["barcode"], records["umi"], genes, records["count"]
            ):
                barcode = int(barcode)
                umi = int(umi)
                gene = int(gene)
                if current_barcode is None:
                    current_barcode, current_umi = barcode, umi
                elif barcode != current_barcode:
                    flush_umi()
                    flush_cell()
                    current_barcode, current_umi = barcode, umi
                elif umi != current_umi:
                    flush_umi()
                    current_umi = umi
                umi_gene_reads[gene] = umi_gene_reads.get(gene, 0) + int(count)
    if current_barcode is not None:
        flush_umi()
        flush_cell()

    if not per_cell_umis or total_reads <= 0:
        raise ValueError(f"No usable reads from filtered cells in {bus_path}")
    umi_matrix = np.vstack(per_cell_umis)
    gene_matrix = np.vstack(per_cell_genes)
    expected_reads = total_reads * fractions
    saturation = 1.0 - np.divide(
        total_expected_umis,
        expected_reads,
        out=np.zeros_like(total_expected_umis),
        where=expected_reads > 0,
    )
    return {
        "cells": len(per_cell_umis),
        "usable_reads": total_reads,
        "expected_reads": expected_reads,
        "expected_umis": total_expected_umis,
        "saturation": saturation,
        "median_umis": np.median(umi_matrix, axis=0),
        "median_genes": np.median(gene_matrix, axis=0),
        "per_cell_umis": umi_matrix,
        "per_cell_genes": gene_matrix,
    }


def build_rows(label, result, fractions):
    return [
        {
            "batch": label,
            "downsample_fraction": float(fraction),
            "cells": int(result["cells"]),
            "expected_usable_reads": float(result["expected_reads"][i]),
            "mean_reads_per_cell": float(result["expected_reads"][i] / result["cells"]),
            "expected_unique_umis": float(result["expected_umis"][i]),
            "sequencing_saturation": float(result["saturation"][i]),
            "median_umis_per_cell": float(result["median_umis"][i]),
            "median_genes_per_cell": float(result["median_genes"][i]),
        }
        for i, fraction in enumerate(fractions)
    ]


def plot_curves(table, output):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for batch, group in table.groupby("batch", sort=False):
        style = {"linewidth": 3, "color": "#1f2937"} if batch == "all" else {"alpha": 0.65}
        x = group["mean_reads_per_cell"]
        axes[0].plot(x, 100 * group["sequencing_saturation"], marker="o", label=batch, **style)
        axes[1].plot(x, group["median_umis_per_cell"], marker="o", **style)
        axes[2].plot(x, group["median_genes_per_cell"], marker="o", **style)
    axes[0].set_ylabel("Sequencing saturation (%)")
    axes[1].set_ylabel("Median UMIs per cell")
    axes[2].set_ylabel("Median genes per cell")
    for ax in axes:
        ax.set_xlabel("Mean usable reads per retained cell")
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("10x-style sequencing-depth rarefaction", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-dirs", nargs="+", required=True)
    parser.add_argument("--filtered-anndata", required=True)
    parser.add_argument("--t2g", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--outdir", default="saturation_qc")
    parser.add_argument(
        "--fractions",
        default="0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
    )
    args = parser.parse_args()
    fractions = np.array([float(x) for x in args.fractions.split(",")], dtype=float)
    if np.any(fractions <= 0) or np.any(fractions > 1) or np.any(np.diff(fractions) <= 0):
        parser.error("--fractions must be strictly increasing values in (0, 1]")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    filtered = ad.read_h5ad(args.filtered_anndata, backed="r")
    covariates = pd.read_csv(args.covariates, dtype=str)
    results = []
    all_results = []
    for directory in sorted(map(Path, args.mapping_dirs), key=lambda path: path.name):
        batch = batch_name(directory)
        bus_path = directory / "output.unfiltered.bus"
        if not bus_path.is_file():
            raise FileNotFoundError(f"Corrected/sorted BUS file not found: {bus_path}")
        ec_to_gene, _gene_count = ec_gene_map(directory, args.t2g)
        result = analyze_bus(bus_path, selected_barcodes(filtered, batch, covariates), ec_to_gene, fractions)
        results.extend(build_rows(batch, result, fractions))
        all_results.append(result)
    filtered.file.close()

    cell_umis = np.vstack([item["per_cell_umis"] for item in all_results])
    cell_genes = np.vstack([item["per_cell_genes"] for item in all_results])
    usable_reads = sum(item["usable_reads"] for item in all_results)
    expected_umis = sum((item["expected_umis"] for item in all_results), np.zeros(len(fractions)))
    expected_reads = usable_reads * fractions
    combined = {
        "cells": cell_umis.shape[0],
        "usable_reads": usable_reads,
        "expected_reads": expected_reads,
        "expected_umis": expected_umis,
        "saturation": 1.0 - expected_umis / expected_reads,
        "median_umis": np.median(cell_umis, axis=0),
        "median_genes": np.median(cell_genes, axis=0),
    }
    results.extend(build_rows("all", combined, fractions))
    table = pd.DataFrame(results)
    table.to_csv(outdir / "sequencing_saturation_curve.tsv", sep="\t", index=False)
    endpoint = table.loc[np.isclose(table["downsample_fraction"], 1.0)].copy()
    endpoint.to_csv(outdir / "sequencing_saturation_metrics.tsv", sep="\t", index=False)
    plot_curves(table, outdir / "sequencing_saturation_curve.png")
    metadata = {
        "method": "10x-style analytic read downsampling from corrected, sorted BUS records",
        "definition": "1 - unique(cell barcode, UMI, gene) molecules / usable reads",
        "cell_population": "cells retained by RNA preprocessing",
        "mapping_filter": "equivalence classes mapping unambiguously to one gene",
        "fractions": fractions.tolist(),
        "mapping_batches": len(all_results),
    }
    (outdir / "sequencing_saturation_method.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
