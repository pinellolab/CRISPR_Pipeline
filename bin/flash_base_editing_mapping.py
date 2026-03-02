#!/usr/bin/env python3

import argparse
import gzip
import json
import os
import shlex
import time
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import torch
except ImportError:  # pragma: no cover - depends on runtime image
    torch = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Notebook-derived base-editing guide mapper with AnnData output compatible with the current pipeline."
    )
    parser.add_argument(
        "-b",
        "--barcode_inclusion_list_fn",
        type=str,
        required=True,
        help="Path to the cell barcode inclusion list file.",
    )
    parser.add_argument(
        "-f",
        "--fastq",
        type=str,
        required=True,
        help="R1/R2 FASTQ files as a whitespace-delimited string: 'R1 R2 [R1 R2 ...]'.",
    )
    parser.add_argument(
        "-g",
        "--guide_set_fn",
        type=str,
        required=True,
        help="Path to the guide set TSV file containing at least 'guide_id' and 'spacer' columns.",
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        type=str,
        default="crispr_map_output",
        help="Output directory prefix.",
    )
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        default=4,
        help="CPU threads to use when running on CPU.",
    )
    parser.add_argument(
        "-d",
        "--downsample_reads",
        type=int,
        default=0,
        help="Limit total read pairs processed across all FASTQ pairs. Set 0 to disable.",
    )
    parser.add_argument(
        "--chemistry",
        type=str,
        required=True,
        help="Parsed seqspec chemistry string in kb format: bc:umi:guide.",
    )
    parser.add_argument(
        "--tolerance",
        type=int,
        default=5,
        help="Maximum allowed Hamming distance for guide matching.",
    )
    parser.add_argument(
        "--guide_len",
        type=int,
        default=0,
        help="Guide length to match. Use 0 to infer from guide metadata.",
    )
    parser.add_argument(
        "--align_on_read",
        type=int,
        default=-1,
        choices=[-1, 0, 1],
        help="Read to align guides against: 0=R1, 1=R2, -1=infer from chemistry.",
    )
    parser.add_argument(
        "--gpu_read_chunk",
        type=int,
        default=4096,
        help="Read chunk size for the dense matcher backend.",
    )
    parser.add_argument(
        "--fastq_chunk_size",
        type=int,
        default=200000,
        help="FASTQ record chunk size for streaming.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Matcher device: 'auto', 'cpu', 'cuda', 'cuda:N', or GPU index like '1'.",
    )
    parser.add_argument(
        "--reverse_complement_guides",
        type=str,
        default="false",
        help="Whether to reverse-complement guide spacers before matching.",
    )
    parser.add_argument(
        "--spacer_tag",
        type=str,
        default="",
        help="Optional prefix added to guide spacers before matching.",
    )
    parser.add_argument(
        "--plot_sample_n",
        type=int,
        default=5000,
        help="Reservoir sample size used for the alignment summary panel.",
    )
    parser.add_argument(
        "--plot_pos_bin_size",
        type=int,
        default=1,
        help="Position bin size for saved alignment plots.",
    )
    parser.add_argument(
        "--plot_read_heatmap_max_rows",
        type=int,
        default=800,
        help="Maximum rows shown in the read-span heatmap.",
    )
    parser.add_argument(
        "--progress_every_chunks",
        type=int,
        default=10,
        help="Emit progress every N FASTQ chunks.",
    )
    return parser.parse_args()


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n", ""}:
        return False
    raise ValueError(f"Cannot parse boolean value from '{value}'.")


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def fastq_seq_iter(path: str) -> Iterator[str]:
    with open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\n")
            handle.readline()
            handle.readline()
            yield seq


def fastq_pair_chunk_iter(
    r1_paths: Sequence[str],
    r2_paths: Sequence[str],
    chunk_size: int,
    max_pairs: Optional[int] = None,
) -> Iterator[Tuple[List[str], List[str]]]:
    remaining = max_pairs
    for r1_path, r2_path in zip(r1_paths, r2_paths):
        it1 = fastq_seq_iter(r1_path)
        it2 = fastq_seq_iter(r2_path)
        while True:
            if remaining is not None and remaining <= 0:
                return
            target = chunk_size if remaining is None else min(chunk_size, remaining)
            r1_seqs, r2_seqs = [], []
            try:
                for _ in range(target):
                    r1_seqs.append(next(it1))
                    r2_seqs.append(next(it2))
            except StopIteration:
                pass
            if not r1_seqs:
                break
            n = min(len(r1_seqs), len(r2_seqs))
            if remaining is not None:
                remaining -= n
            yield r1_seqs[:n], r2_seqs[:n]
            if n < target:
                break


def parse_fastq_pairs(fastq_arg: str) -> Tuple[List[str], List[str]]:
    tokens = shlex.split(fastq_arg)
    if len(tokens) < 2 or len(tokens) % 2 != 0:
        raise ValueError(
            "FASTQ input must contain an even number of paths formatted as 'R1 R2 [R1 R2 ...]'."
        )
    return tokens[0::2], tokens[1::2]


def parse_triplet(spec: str) -> Tuple[int, int, int]:
    fields = [x.strip() for x in spec.split(",")]
    if len(fields) != 3:
        raise ValueError(f"Invalid chemistry triplet '{spec}'. Expected 'read,start,end'.")
    return int(fields[0]), int(fields[1]), int(fields[2])


def parse_seqspec_chemistry(chemistry: str) -> Tuple[Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]]:
    parts = [p.strip() for p in chemistry.strip().split(":") if p.strip()]
    if len(parts) != 3:
        raise ValueError(
            f"Invalid chemistry string '{chemistry}'. Expected 3 colon-separated triplets: bc:umi:guide."
        )
    return parse_triplet(parts[0]), parse_triplet(parts[1]), parse_triplet(parts[2])


def slice_from_read_pair(seq_r1: str, seq_r2: str, spec: Tuple[int, int, int]) -> str:
    read_idx, start, end = spec
    src = seq_r1 if read_idx == 0 else seq_r2
    return src[start:end]


_DNA2BIT = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
_BASE2INT = np.full(256, 4, dtype=np.uint8)
for base, value in [(b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3), (b"a", 0), (b"c", 1), (b"g", 2), (b"t", 3)]:
    _BASE2INT[base[0]] = value


def pack_dna_2bit(seq: str) -> Optional[int]:
    packed = 0
    for char in seq:
        value = _DNA2BIT.get(char)
        if value is None:
            return None
        packed = (packed << 2) | value
    return packed


def encode_seqs_uint8(seqs: Sequence[str]) -> List[np.ndarray]:
    encoded = []
    for seq in seqs:
        seq_bytes = str(seq).encode("ascii", "ignore")
        encoded.append(_BASE2INT[np.frombuffer(seq_bytes, dtype=np.uint8)])
    return encoded


def pad_reads_u8(encoded_list: Sequence[np.ndarray], pad_value: int = 4):
    lengths = np.array([len(x) for x in encoded_list], dtype=np.int32)
    max_len = int(lengths.max()) if len(lengths) else 0
    matrix = np.full((len(encoded_list), max_len), pad_value, dtype=np.uint8)
    for idx, arr in enumerate(encoded_list):
        matrix[idx, : len(arr)] = arr
    return matrix, lengths, max_len


def resolve_device(requested: str) -> str:
    text = str(requested).strip()
    lowered = text.lower()
    if lowered == "auto":
        if torch is not None and torch.cuda.is_available():
            return "cuda"
        return "cpu"
    if lowered.isdigit():
        if torch is None or not torch.cuda.is_available():
            raise RuntimeError(f"Device '{requested}' requested but CUDA is not available.")
        return f"cuda:{lowered}"
    if lowered.startswith("cuda"):
        if torch is None or not torch.cuda.is_available():
            raise RuntimeError(f"Device '{requested}' requested but CUDA is not available.")
        return lowered
    return "cpu"


def gpu_mem_str(device: str) -> str:
    if torch is None or not device.startswith("cuda"):
        return ""
    dev = torch.device(device)
    alloc = torch.cuda.memory_allocated(dev) / (1024**3)
    reserved = torch.cuda.memory_reserved(dev) / (1024**3)
    free, total = torch.cuda.mem_get_info(dev)
    used = (total - free) / (1024**3)
    return f" | GPU mem GB: alloc={alloc:.2f} reserved={reserved:.2f} driver_used={used:.2f}/{total / (1024**3):.2f}"


class DenseGuideMatcher:
    def __init__(self, guides: Sequence[str], guide_len: int, tolerance: int, device: str, read_chunk: int):
        self.guides = list(guides)
        self.L = int(guide_len)
        self.tolerance = int(tolerance)
        self.device = resolve_device(device)
        self.read_chunk = int(read_chunk)
        self.guide_mat = np.stack([arr[: self.L] for arr in encode_seqs_uint8(self.guides)], axis=0).astype(np.uint8)
        self.torch_backend = torch is not None

        if self.torch_backend:
            self.dev = torch.device(self.device)
            self.guide_t = torch.from_numpy(self.guide_mat).to(device=self.dev, dtype=torch.uint8)

    def best_hits(self, reads: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.torch_backend:
            return self._best_hits_torch(reads)
        return self._best_hits_numpy(reads)

    def _best_hits_torch(self, reads: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        out_dist, out_start, out_gid = [], [], []
        n_reads = len(reads)
        guide_count = self.guide_t.shape[0]

        for start_idx in range(0, n_reads, self.read_chunk):
            chunk = reads[start_idx : start_idx + self.read_chunk]
            encoded_reads = encode_seqs_uint8(chunk)
            read_mat, lengths, max_len = pad_reads_u8(encoded_reads, pad_value=4)
            batch_size = read_mat.shape[0]

            if max_len < self.L:
                out_dist.append(np.full((batch_size,), 10**9, dtype=np.int32))
                out_start.append(np.full((batch_size,), -1, dtype=np.int32))
                out_gid.append(np.full((batch_size,), -1, dtype=np.int32))
                continue

            read_t = torch.from_numpy(read_mat).to(device=self.dev, dtype=torch.uint8)
            lens_t = torch.from_numpy(lengths).to(device=self.dev, dtype=torch.int32)
            span_count = max_len - self.L + 1

            windows = read_t.unfold(1, self.L, 1)
            max_start_each = (lens_t - self.L).clamp(min=-1)
            starts = torch.arange(0, span_count, device=self.dev, dtype=torch.int32)[None, :]
            valid = starts <= max_start_each[:, None]

            mismatches = torch.sum(
                windows[:, :, None, :] != self.guide_t[None, None, :, :],
                dim=-1,
                dtype=torch.int16,
            )
            mismatches = torch.where(
                valid[:, :, None],
                mismatches,
                torch.tensor(32767, device=self.dev, dtype=torch.int16),
            )

            flattened = mismatches.reshape(batch_size, span_count * guide_count)
            best_dist_i16, argmin = torch.min(flattened, dim=1)
            best_start_i32 = (argmin // guide_count).to(torch.int32)
            best_gid_i32 = (argmin % guide_count).to(torch.int32)

            ok = best_dist_i16 <= self.tolerance
            best_dist_i32 = best_dist_i16.to(torch.int32)
            best_dist_i32 = torch.where(ok, best_dist_i32, torch.tensor(10**9, device=self.dev, dtype=torch.int32))
            best_start_i32 = torch.where(ok, best_start_i32, torch.tensor(-1, device=self.dev, dtype=torch.int32))
            best_gid_i32 = torch.where(ok, best_gid_i32, torch.tensor(-1, device=self.dev, dtype=torch.int32))

            out_dist.append(best_dist_i32.cpu().numpy())
            out_start.append(best_start_i32.cpu().numpy())
            out_gid.append(best_gid_i32.cpu().numpy())

        return np.concatenate(out_dist), np.concatenate(out_start), np.concatenate(out_gid)

    def _best_hits_numpy(self, reads: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        out_dist, out_start, out_gid = [], [], []
        n_reads = len(reads)

        for start_idx in range(0, n_reads, self.read_chunk):
            chunk = reads[start_idx : start_idx + self.read_chunk]
            encoded_reads = encode_seqs_uint8(chunk)
            read_mat, lengths, max_len = pad_reads_u8(encoded_reads, pad_value=4)
            batch_size = read_mat.shape[0]

            if max_len < self.L:
                out_dist.append(np.full((batch_size,), 10**9, dtype=np.int32))
                out_start.append(np.full((batch_size,), -1, dtype=np.int32))
                out_gid.append(np.full((batch_size,), -1, dtype=np.int32))
                continue

            span_count = max_len - self.L + 1
            windows = np.lib.stride_tricks.sliding_window_view(read_mat, self.L, axis=1)
            valid = np.arange(span_count, dtype=np.int32)[None, :] <= (lengths - self.L)[:, None]

            best_dist = np.full(batch_size, 10**9, dtype=np.int32)
            best_start = np.full(batch_size, -1, dtype=np.int32)
            best_gid = np.full(batch_size, -1, dtype=np.int32)

            for gid, guide in enumerate(self.guide_mat):
                mism = np.sum(windows != guide[None, None, :], axis=-1, dtype=np.int16)
                mism = np.where(valid, mism, np.int16(32767))
                local_start = mism.argmin(axis=1).astype(np.int32)
                local_dist = mism[np.arange(batch_size), local_start].astype(np.int32)

                better = local_dist < best_dist
                better_start = (local_dist == best_dist) & (
                    (best_start < 0) | (local_start < best_start)
                )
                better_gid = (local_dist == best_dist) & (local_start == best_start) & (
                    (best_gid < 0) | (gid < best_gid)
                )
                update = better | better_start | better_gid

                best_dist[update] = local_dist[update]
                best_start[update] = local_start[update]
                best_gid[update] = gid

            ok = best_dist <= self.tolerance
            best_dist = np.where(ok, best_dist, 10**9).astype(np.int32)
            best_start = np.where(ok, best_start, -1).astype(np.int32)
            best_gid = np.where(ok, best_gid, -1).astype(np.int32)

            out_dist.append(best_dist)
            out_start.append(best_start)
            out_gid.append(best_gid)

        return np.concatenate(out_dist), np.concatenate(out_start), np.concatenate(out_gid)


def add_umi(cell_map: Dict[int, Dict[int, set]], cell_key: int, guide_key: int, umi_key: int):
    guide_map = cell_map.get(cell_key)
    if guide_map is None:
        guide_map = {}
        cell_map[cell_key] = guide_map
    umi_set = guide_map.get(guide_key)
    if umi_set is None:
        umi_set = set()
        guide_map[guide_key] = umi_set
    umi_set.add(umi_key)


class TimerLog:
    def __init__(self):
        self.records = []

    @contextmanager
    def section(self, name: str):
        start = time.time()
        yield
        self.records.append((name, time.time() - start))

    def summary(self) -> pd.DataFrame:
        df = pd.DataFrame(self.records, columns=["step", "seconds"])
        if df.empty:
            return df
        df = df.sort_values("seconds", ascending=False).reset_index(drop=True)
        total = df["seconds"].sum()
        df["pct"] = df["seconds"] / total * 100 if total > 0 else 0.0
        return df


@dataclass
class PlotStats:
    tol: int
    guide_len: int
    n_guides: int
    pos_bin_size: int = 1
    sample_n: int = 5000
    sample_seed: int = 0
    total_reads: int = 0
    no_hit: int = 0
    dist_counts: Optional[np.ndarray] = None
    mismatch_pos_counts: Optional[np.ndarray] = None
    read_base_at_mismatch: Optional[np.ndarray] = None
    guide_base_at_mismatch: Optional[np.ndarray] = None
    sub_matrix: Optional[np.ndarray] = None
    guide_counts: Optional[np.ndarray] = None
    rng: Optional[np.random.Generator] = None
    sample_filled: int = 0
    sample_seen: int = 0
    sample_read_len: Optional[np.ndarray] = None
    sample_start: Optional[np.ndarray] = None
    sample_end: Optional[np.ndarray] = None
    sample_gid: Optional[np.ndarray] = None
    max_read_len_seen: int = 0

    def __post_init__(self):
        self.dist_counts = np.zeros(self.tol + 1, dtype=np.int64)
        self.mismatch_pos_counts = np.zeros(self.guide_len, dtype=np.int64)
        self.read_base_at_mismatch = np.zeros((4, self.guide_len), dtype=np.int64)
        self.guide_base_at_mismatch = np.zeros((4, self.guide_len), dtype=np.int64)
        self.sub_matrix = np.zeros((4, 4), dtype=np.int64)
        self.guide_counts = np.zeros(self.n_guides, dtype=np.int64)
        self.rng = np.random.default_rng(self.sample_seed)
        self.sample_read_len = np.zeros(self.sample_n, dtype=np.int32)
        self.sample_start = np.full(self.sample_n, -1, dtype=np.int32)
        self.sample_end = np.full(self.sample_n, -1, dtype=np.int32)
        self.sample_gid = np.full(self.sample_n, -1, dtype=np.int32)

    @staticmethod
    def _base_to_idx(ch: str) -> int:
        if ch == "A":
            return 0
        if ch == "C":
            return 1
        if ch == "G":
            return 2
        if ch == "T":
            return 3
        return -1

    def _maybe_store_sample(self, read_len: int, start: int, end: int, gid: int):
        self.sample_seen += 1
        if self.sample_filled < self.sample_n:
            idx = self.sample_filled
            self.sample_filled += 1
        else:
            replacement = self.rng.integers(0, self.sample_seen)
            if replacement >= self.sample_n:
                return
            idx = int(replacement)

        self.sample_read_len[idx] = read_len
        self.sample_start[idx] = start
        self.sample_end[idx] = end
        self.sample_gid[idx] = gid

    def update(self, reads: Sequence[str], guides: Sequence[str], best_dist, best_start, best_gid):
        self.total_reads += len(reads)
        self.max_read_len_seen = max(self.max_read_len_seen, max((len(r) for r in reads), default=0))

        for read, dist, start, gid in zip(reads, best_dist, best_start, best_gid):
            if gid < 0 or start < 0 or dist >= 10**8:
                self.no_hit += 1
                self._maybe_store_sample(len(read), -1, -1, -1)
                continue

            dist = int(dist)
            start = int(start)
            gid = int(gid)
            end = start + self.guide_len
            if dist <= self.tol:
                self.dist_counts[dist] += 1
            else:
                self.no_hit += 1
                self._maybe_store_sample(len(read), -1, -1, -1)
                continue

            self.guide_counts[gid] += 1
            guide_seq = guides[gid]
            read_window = read[start:end]

            for pos, (read_base, guide_base) in enumerate(zip(read_window, guide_seq)):
                if read_base == guide_base:
                    continue
                self.mismatch_pos_counts[pos] += 1
                read_idx = self._base_to_idx(read_base)
                guide_idx = self._base_to_idx(guide_base)
                if read_idx >= 0:
                    self.read_base_at_mismatch[read_idx, pos] += 1
                if guide_idx >= 0:
                    self.guide_base_at_mismatch[guide_idx, pos] += 1
                if read_idx >= 0 and guide_idx >= 0:
                    self.sub_matrix[guide_idx, read_idx] += 1

            self._maybe_store_sample(len(read), start, end, gid)


def downsample_rows_mean(matrix: np.ndarray, target_rows: int):
    if matrix.shape[0] <= target_rows or target_rows <= 0:
        return matrix
    edges = np.linspace(0, matrix.shape[0], target_rows + 1).astype(int)
    out = np.zeros((target_rows, matrix.shape[1]), dtype=float)
    for idx in range(target_rows):
        left, right = edges[idx], edges[idx + 1]
        if right <= left:
            right = min(matrix.shape[0], left + 1)
        out[idx] = matrix[left:right].mean(axis=0)
    return out


def save_alignment_panel(
    stats: PlotStats,
    guides: Sequence[str],
    output_path: str,
    pos_bin_size: int = 1,
    read_heatmap_max_rows: int = 800,
):
    total = stats.total_reads
    if total == 0:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No reads processed", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    exact = int(stats.dist_counts[0]) if stats.dist_counts is not None and len(stats.dist_counts) else 0
    no_hit = int(stats.no_hit)
    mismatched_any = max(0, total - exact - no_hit)
    used_hits = max(0, total - no_hit)
    guide_counts_s = pd.Series(stats.guide_counts, index=np.arange(stats.n_guides)).sort_values(ascending=False)

    sample_count = stats.sample_filled
    sample_start = stats.sample_start[:sample_count]
    sample_end = stats.sample_end[:sample_count]
    sample_gid = stats.sample_gid[:sample_count]

    max_read_len = max(1, stats.max_read_len_seen)
    max_start = max(0, max_read_len - stats.guide_len)
    n_bins = (max_start // pos_bin_size) + 1
    guide_start_heatmap = np.zeros((stats.n_guides, n_bins), dtype=np.int32)

    valid = (sample_gid >= 0) & (sample_start >= 0)
    if np.any(valid):
        gids = sample_gid[valid].astype(np.int32, copy=False)
        starts = np.clip((sample_start[valid] // pos_bin_size).astype(np.int32, copy=False), 0, n_bins - 1)
        np.add.at(guide_start_heatmap, (gids, starts), 1)

    n_pos_bins = (max_read_len // pos_bin_size) + 1
    order = np.argsort(np.where(valid, sample_start, 10**12))
    read_span_heatmap = np.zeros((sample_count, n_pos_bins), dtype=float)
    for row_idx, sample_idx in enumerate(order):
        if sample_start[sample_idx] < 0 or sample_end[sample_idx] < 0:
            continue
        left = max(0, int(sample_start[sample_idx] // pos_bin_size))
        right = min(n_pos_bins, int((sample_end[sample_idx] - 1) // pos_bin_size + 1))
        if right > left:
            read_span_heatmap[row_idx, left:right] = 1.0

    read_span_plot = downsample_rows_mean(read_span_heatmap, read_heatmap_max_rows)
    bases = ["A", "C", "G", "T"]
    step = max(1, stats.guide_len // 10)

    fig, axes = plt.subplots(3, 4, figsize=(28, 14))

    pie_sizes = [exact, mismatched_any, no_hit]
    pie_labels = [
        f"Exact (0)\n{exact}/{total}",
        f"Mismatched (1-{stats.tol})\n{mismatched_any}/{total}",
        f"No hit\n{no_hit}/{total}",
    ]
    axes[0, 0].pie(pie_sizes, labels=pie_labels, autopct=lambda p: f"{p:.1f}%")
    axes[0, 0].set_title(f"Read outcomes (tol={stats.tol})")

    x_labels = [str(k) for k in range(stats.tol + 1)] + ["no_hit"]
    y_vals = [int(stats.dist_counts[k]) for k in range(stats.tol + 1)] + [no_hit]
    x_vals = np.arange(len(x_labels))
    axes[0, 1].bar(x_vals, y_vals)
    axes[0, 1].set_xticks(x_vals)
    axes[0, 1].set_xticklabels(x_labels)
    axes[0, 1].set_xlabel("Best mismatches per read")
    axes[0, 1].set_ylabel("Reads")
    axes[0, 1].set_title("Counts by mismatch number + no hit")

    pos_x = np.arange(1, stats.guide_len + 1)
    axes[0, 2].bar(pos_x, stats.mismatch_pos_counts)
    axes[0, 2].set_xlabel("Position in aligned guide (1-based)")
    axes[0, 2].set_ylabel("Mismatch count")
    axes[0, 2].set_title(f"Mismatches per position (n_hits={used_hits})")
    axes[0, 2].set_xticks(np.arange(1, stats.guide_len + 1, step))

    if guide_counts_s.sum() == 0:
        axes[0, 3].text(0.5, 0.5, "No guide hits", ha="center", va="center")
        axes[0, 3].set_axis_off()
    else:
        axes[0, 3].bar(np.arange(len(guide_counts_s)), guide_counts_s.values)
        axes[0, 3].set_xlabel("Guide (sorted by count)")
        axes[0, 3].set_ylabel("Read count")
        axes[0, 3].set_title("All guides by count")
        axes[0, 3].set_xticks([])

    im = axes[1, 0].imshow(stats.read_base_at_mismatch, aspect="auto")
    axes[1, 0].set_yticks([0, 1, 2, 3])
    axes[1, 0].set_yticklabels(bases)
    axes[1, 0].set_xticks(np.arange(0, stats.guide_len, step))
    axes[1, 0].set_xticklabels([str(v) for v in range(1, stats.guide_len + 1, step)])
    axes[1, 0].set_title("Read base at mismatches")
    fig.colorbar(im, ax=axes[1, 0], fraction=0.046, pad=0.04)

    im = axes[1, 1].imshow(stats.guide_base_at_mismatch, aspect="auto")
    axes[1, 1].set_yticks([0, 1, 2, 3])
    axes[1, 1].set_yticklabels(bases)
    axes[1, 1].set_xticks(np.arange(0, stats.guide_len, step))
    axes[1, 1].set_xticklabels([str(v) for v in range(1, stats.guide_len + 1, step)])
    axes[1, 1].set_title("Guide base at mismatches")
    fig.colorbar(im, ax=axes[1, 1], fraction=0.046, pad=0.04)

    im = axes[1, 2].imshow(stats.sub_matrix, aspect="auto")
    axes[1, 2].set_xticks([0, 1, 2, 3])
    axes[1, 2].set_xticklabels(bases)
    axes[1, 2].set_yticks([0, 1, 2, 3])
    axes[1, 2].set_yticklabels(bases)
    axes[1, 2].set_xlabel("Read base")
    axes[1, 2].set_ylabel("Guide base")
    axes[1, 2].set_title("Substitution matrix")
    fig.colorbar(im, ax=axes[1, 2], fraction=0.046, pad=0.04)

    if guide_counts_s.sum() == 0:
        axes[1, 3].text(0.5, 0.5, "No guide hits", ha="center", va="center")
        axes[1, 3].set_axis_off()
    else:
        ranks = np.arange(1, len(guide_counts_s) + 1)
        axes[1, 3].plot(ranks, guide_counts_s.values, marker="o", linestyle="-")
        axes[1, 3].set_xlabel("Guide rank")
        axes[1, 3].set_ylabel("Read count")
        axes[1, 3].set_title("Guide rank plot")

    if guide_start_heatmap.sum() == 0:
        axes[2, 0].text(0.5, 0.5, "No alignments to plot", ha="center", va="center")
        axes[2, 0].set_axis_off()
    else:
        guide_order = guide_counts_s.index.to_numpy(dtype=int)
        im = axes[2, 0].imshow(guide_start_heatmap[guide_order, :], aspect="auto")
        axes[2, 0].set_ylabel("Guide (sorted by count)")
        axes[2, 0].set_xlabel(f"Start position bin (bin={pos_bin_size} bp)")
        axes[2, 0].set_title(f"Guide alignment start positions (sampled {sample_count} reads)")
        axes[2, 0].set_yticks([])
        fig.colorbar(im, ax=axes[2, 0], fraction=0.046, pad=0.04)

    if np.nansum(read_span_plot) == 0:
        axes[2, 1].text(0.5, 0.5, "No read spans to plot", ha="center", va="center")
        axes[2, 1].set_axis_off()
    else:
        im = axes[2, 1].imshow(read_span_plot, aspect="auto")
        axes[2, 1].set_xlabel(f"Read position bin (bin={pos_bin_size} bp)")
        axes[2, 1].set_ylabel("Reads")
        axes[2, 1].set_title("Aligned guide span per read")
        fig.colorbar(im, ax=axes[2, 1], fraction=0.046, pad=0.04)

    axes[2, 2].set_axis_off()
    axes[2, 3].set_axis_off()

    fig.suptitle(
        f"Alignment summary (tol={stats.tol}; L={stats.guide_len}; n_reads={total}; sampled={sample_count})",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def load_guide_metadata(path: str, reverse_complement_guides: bool, spacer_tag: str):
    guide_df = pd.read_csv(path, sep="\t")
    required = {"guide_id", "spacer"}
    missing = required.difference(guide_df.columns)
    if missing:
        raise ValueError(f"Guide metadata is missing required columns: {sorted(missing)}")

    raw_spacers = guide_df["spacer"].astype(str).tolist()
    match_spacers = raw_spacers.copy()
    if reverse_complement_guides:
        match_spacers = [reverse_complement(seq) for seq in match_spacers]
    if spacer_tag:
        match_spacers = [f"{spacer_tag}{seq}" for seq in match_spacers]

    return guide_df, raw_spacers, match_spacers


def determine_guide_len(match_spacers: Sequence[str], requested_len: int) -> int:
    if requested_len and requested_len > 0:
        too_short = [seq for seq in match_spacers if len(seq) < requested_len]
        if too_short:
            raise ValueError("Requested guide_len is longer than at least one guide sequence.")
        return requested_len

    lengths = sorted({len(seq) for seq in match_spacers})
    if len(lengths) != 1:
        raise ValueError(
            f"Guide sequences have multiple lengths {lengths}. Set --guide_len explicitly for flash mapping."
        )
    return lengths[0]


def build_unique_guide_map(match_spacers: Sequence[str], guide_len: int):
    unique_guides = []
    seen = {}
    output_to_unique = []
    for seq in match_spacers:
        clipped = seq[:guide_len]
        if clipped not in seen:
            seen[clipped] = len(unique_guides)
            unique_guides.append(clipped)
        output_to_unique.append(seen[clipped])
    return unique_guides, np.asarray(output_to_unique, dtype=np.int32)


def load_barcode_whitelist(path: str) -> Tuple[set, int]:
    whitelist = set()
    total = 0
    with open_maybe_gzip(path) as handle:
        for line in handle:
            barcode = line.strip().split("\t")[0]
            if not barcode:
                continue
            total += 1
            packed = pack_dna_2bit(barcode)
            if packed is not None:
                whitelist.add(packed)
    return whitelist, total


def build_output_anndata(
    counts_by_cell: Dict[int, Dict[int, set]],
    barcode_key_to_str: Dict[int, str],
    guide_df: pd.DataFrame,
    raw_spacers: Sequence[str],
    output_to_unique: np.ndarray,
):
    cell_keys = np.array(sorted(counts_by_cell.keys()), dtype=np.int64)
    n_cells = len(cell_keys)
    n_guides_unique = int(output_to_unique.max()) + 1 if len(output_to_unique) else 0

    rows, cols, vals = [], [], []
    cell_index = {cell_key: idx for idx, cell_key in enumerate(cell_keys)}
    for cell_key, guide_map in counts_by_cell.items():
        row_idx = cell_index[cell_key]
        for unique_gid, umi_set in guide_map.items():
            rows.append(row_idx)
            cols.append(unique_gid)
            vals.append(len(umi_set))

    x_unique = sp.csr_matrix((vals, (rows, cols)), shape=(n_cells, n_guides_unique), dtype=np.int32)
    x_output = x_unique[:, output_to_unique] if len(output_to_unique) else sp.csr_matrix((n_cells, 0), dtype=np.int32)

    adata = ad.AnnData(
        X=x_output,
        obs=pd.DataFrame(index=[barcode_key_to_str[int(cell_key)] for cell_key in cell_keys]),
        var=pd.DataFrame(index=pd.Index([str(seq) for seq in raw_spacers])),
    )
    adata.var["guide_id"] = guide_df["guide_id"].astype(str).values
    return adata


def write_dashboard_json(
    output_dir: str,
    n_pairs: int,
    n_assigned: int,
    total_umis: int,
    adata: ad.AnnData,
    onlist_read_pairs: int,
    onlist_barcodes: int,
    whitelist_size: int,
    guide_count_output: int,
):
    n_cells = int(adata.n_obs)
    processed = max(1, int(n_pairs))
    inspect = {
        "numRecords": int(n_pairs),
        "numReads": int(n_pairs),
        "numBarcodes": n_cells,
        "numUMIs": int(total_umis),
        "numBarcodeUMIs": int(adata.X.nnz),
        "gtRecords": int(guide_count_output),
        "numBarcodesOnOnlist": int(onlist_barcodes),
        "numReadsOnOnlist": int(onlist_read_pairs),
        "percentageBarcodesOnOnlist": (100.0 * onlist_barcodes / max(1, n_cells)),
        "percentageReadsOnOnlist": (100.0 * onlist_read_pairs / processed),
        "meanReadsPerBarcode": float(n_assigned / max(1, n_cells)),
        "meanUMIsPerBarcode": float(total_umis / max(1, n_cells)),
        "barcodeWhitelistSize": int(whitelist_size),
    }
    run_info = {
        "n_targets": int(guide_count_output),
        "n_processed": int(n_pairs),
        "n_unique": int(total_umis),
        "n_pseudoaligned": int(n_assigned),
        "p_pseudoaligned": (100.0 * n_assigned / processed),
        "p_unique": (100.0 * total_umis / processed),
        "call": "flash_base_editing_mapping.py",
    }

    with open(os.path.join(output_dir, "inspect.json"), "w") as handle:
        json.dump(inspect, handle, indent=2)
    with open(os.path.join(output_dir, "run_info.json"), "w") as handle:
        json.dump(run_info, handle, indent=2)


def main():
    args = parse_args()
    reverse_complement_guides = parse_bool(args.reverse_complement_guides)
    r1_paths, r2_paths = parse_fastq_pairs(args.fastq)
    barcode_spec, umi_spec, guide_spec = parse_seqspec_chemistry(args.chemistry)
    align_on_read = guide_spec[0] if args.align_on_read < 0 else args.align_on_read
    downsample_reads = None if args.downsample_reads <= 0 else args.downsample_reads

    if torch is not None and resolve_device(args.device) == "cpu":
        torch.set_num_threads(max(1, int(args.cores)))

    timer = TimerLog()
    output_counts_dir = os.path.join(args.output_prefix, "counts_unfiltered")
    output_plot_dir = os.path.join(args.output_prefix, "plots")
    os.makedirs(output_counts_dir, exist_ok=True)
    os.makedirs(output_plot_dir, exist_ok=True)

    with timer.section("load inputs"):
        guide_df, raw_spacers, match_spacers = load_guide_metadata(
            args.guide_set_fn,
            reverse_complement_guides=reverse_complement_guides,
            spacer_tag=args.spacer_tag,
        )
        guide_len = determine_guide_len(match_spacers, args.guide_len)
        unique_guides, output_to_unique = build_unique_guide_map(match_spacers, guide_len)
        barcode_whitelist, whitelist_size = load_barcode_whitelist(args.barcode_inclusion_list_fn)
        matcher = DenseGuideMatcher(
            guides=unique_guides,
            guide_len=guide_len,
            tolerance=args.tolerance,
            device=args.device,
            read_chunk=args.gpu_read_chunk,
        )
        stats = PlotStats(
            tol=args.tolerance,
            guide_len=guide_len,
            n_guides=len(unique_guides),
            pos_bin_size=args.plot_pos_bin_size,
            sample_n=args.plot_sample_n,
        )

    barcode_key_to_str: Dict[int, str] = {}
    counts_by_cell: Dict[int, Dict[int, set]] = {}
    n_pairs = 0
    n_assigned = 0
    n_nohit = 0
    n_bad_bc = 0
    n_bad_umi = 0
    onlist_read_pairs = 0

    print(
        f"Starting flash base-editing mapping with {len(unique_guides)} unique guide sequences "
        f"expanded to {len(raw_spacers)} output guides. device={matcher.device}"
    )

    start_time = time.time()
    for chunk_idx, (r1_seqs, r2_seqs) in enumerate(
        fastq_pair_chunk_iter(r1_paths, r2_paths, args.fastq_chunk_size, downsample_reads),
        start=1,
    ):
        with timer.section(f"chunk_{chunk_idx:05d}"):
            chunk_size = len(r1_seqs)
            n_pairs += chunk_size

            bc_keys = np.empty(chunk_size, dtype=np.int64)
            umi_keys = np.empty(chunk_size, dtype=np.int64)
            valid_bu = np.ones(chunk_size, dtype=bool)
            onlist_mask = np.zeros(chunk_size, dtype=bool)

            for idx, (seq_r1, seq_r2) in enumerate(zip(r1_seqs, r2_seqs)):
                barcode = slice_from_read_pair(seq_r1, seq_r2, barcode_spec)
                umi = slice_from_read_pair(seq_r1, seq_r2, umi_spec)

                barcode_key = pack_dna_2bit(barcode)
                umi_key = pack_dna_2bit(umi)
                if barcode_key is None:
                    valid_bu[idx] = False
                    n_bad_bc += 1
                    continue
                if umi_key is None:
                    valid_bu[idx] = False
                    n_bad_umi += 1
                    continue

                bc_keys[idx] = barcode_key
                umi_keys[idx] = umi_key
                barcode_key_to_str.setdefault(barcode_key, barcode)
                if not barcode_whitelist or barcode_key in barcode_whitelist:
                    onlist_mask[idx] = True
                    onlist_read_pairs += 1

            align_reads = r2_seqs if align_on_read == 1 else r1_seqs
            best_dist, best_start, best_gid = matcher.best_hits(align_reads)
            stats.update(align_reads, unique_guides, best_dist, best_start, best_gid)

            for idx in range(chunk_size):
                if not valid_bu[idx]:
                    continue
                gid = int(best_gid[idx])
                if gid < 0:
                    n_nohit += 1
                    continue
                n_assigned += 1
                add_umi(counts_by_cell, int(bc_keys[idx]), gid, int(umi_keys[idx]))

        if chunk_idx % max(1, args.progress_every_chunks) == 0:
            elapsed = time.time() - start_time
            rate = n_pairs / max(elapsed, 1e-9)
            print(
                f"[chunk {chunk_idx}] pairs={n_pairs:,} assigned={n_assigned:,} no_hit={n_nohit:,} "
                f"bad_bc={n_bad_bc:,} bad_umi={n_bad_umi:,} | {rate:,.0f} pairs/s{gpu_mem_str(matcher.device)}"
            )

    with timer.section("build outputs"):
        adata = build_output_anndata(
            counts_by_cell=counts_by_cell,
            barcode_key_to_str=barcode_key_to_str,
            guide_df=guide_df,
            raw_spacers=raw_spacers,
            output_to_unique=output_to_unique,
        )
        adata.uns["params"] = {
            "fastq_pairs": list(zip(r1_paths, r2_paths)),
            "chemistry": args.chemistry,
            "guide_len": guide_len,
            "tolerance": args.tolerance,
            "align_on_read": align_on_read,
            "fastq_chunk_size": args.fastq_chunk_size,
            "gpu_read_chunk": args.gpu_read_chunk,
            "requested_device": args.device,
            "resolved_device": matcher.device,
            "reverse_complement_guides": reverse_complement_guides,
            "spacer_tag": args.spacer_tag,
        }
        adata.uns["qc"] = {
            "total_pairs": int(n_pairs),
            "assigned": int(n_assigned),
            "no_hit": int(n_nohit),
            "bad_bc": int(n_bad_bc),
            "bad_umi": int(n_bad_umi),
            "n_cells": int(adata.n_obs),
            "n_guides_output": int(adata.n_vars),
            "n_guides_aligned_unique": int(len(unique_guides)),
        }
        output_h5ad = os.path.join(output_counts_dir, "adata.h5ad")
        adata.write_h5ad(output_h5ad)

        total_umis = int(adata.X.sum())
        onlist_barcodes = int(
            sum(1 for barcode in adata.obs_names if pack_dna_2bit(barcode) in barcode_whitelist)
        ) if barcode_whitelist else int(adata.n_obs)

        write_dashboard_json(
            output_dir=args.output_prefix,
            n_pairs=n_pairs,
            n_assigned=n_assigned,
            total_umis=total_umis,
            adata=adata,
            onlist_read_pairs=onlist_read_pairs if barcode_whitelist else n_pairs,
            onlist_barcodes=onlist_barcodes,
            whitelist_size=whitelist_size,
            guide_count_output=adata.n_vars,
        )

    with timer.section("plots"):
        save_alignment_panel(
            stats=stats,
            guides=unique_guides,
            output_path=os.path.join(output_plot_dir, "alignment_summary_panel.png"),
            pos_bin_size=args.plot_pos_bin_size,
            read_heatmap_max_rows=args.plot_read_heatmap_max_rows,
        )

    timings = timer.summary()
    timings.to_csv(os.path.join(args.output_prefix, "timings.tsv"), sep="\t", index=False)

    print("\nTiming summary:")
    print(timings)
    print(f"\nSaved AnnData: {os.path.join(output_counts_dir, 'adata.h5ad')}")
    print(adata)


if __name__ == "__main__":
    main()
