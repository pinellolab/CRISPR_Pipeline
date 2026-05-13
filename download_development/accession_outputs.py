#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Register the standard first five Perturb-seq pipeline outputs (matrix + four tabular
files) on the IGVF portal by invoking ``iu_register`` automatically (temporary TSV
inputs only; no artifact written to disk beyond logs).

Schemas follow igvfd (matrix_file, tabular_file). Content types match portal enums:
  - inference *.h5mu -> matrix_file, content_type "annotated multimodal CRISPR matrix"
  - *per_guide*.tsv.gz -> tabular_file, "differential guide quantifications"
  - *per_element*.tsv.gz -> tabular_file, "differential element quantifications"

reference_files are resolved from the pipeline params and portal metadata:
  - if ``pipeline_info/params*.json`` says ``use_igvf_reference=true`` (default),
    the script uses IGVF reference ``IGVFFI9561BASO`` and submits its
    own accession plus ``derived_from`` values as ``reference_files``
  - if ``use_igvf_reference`` is false, pass ``--reference-files`` explicitly
    (comma-separated accessions or aliases)

ReferenceFile lookups use the same API host as submissions: ``-m prod`` →
``api.data.igvf.org``, ``-m staging`` → ``api.staging.igvf.org`` (or ``IGVF_MODE``
when ``-m`` is omitted). Override only if needed with ``--portal-api``.

**Re-upload mode:** add ``--reupload`` to push GCS objects onto
existing File records on an analysis set (md5-matched; errors if the same md5 exists on
another file set, or if the set contains a File whose md5 is not one of the pipeline
outputs, or a File with no md5sum). Requires igvf_utils ``Connection``, API keys, and
AWS keys per igvf_utils.

Requires: ``igvf_utils`` with CLI/registrar for POST. If input is ``gs://``, requires
``gsutil`` on PATH for listing/params read; for md5sum on gs:// objects without hashing
via registrar deps, pass ``--use-gsutil-hash``.

After building rows, syncs with the portal: for each pipeline md5 that already exists on
the analysis set, if ``upload_status`` is neither ``validated`` nor ``validation exempted``,
the GCS object is uploaded automatically via ``Connection.upload_file``. Remaining files (no matching md5 on the file set)
are POSTed via ``iu_register``. If the set contains Files whose md5 does not match any of
the five outputs, or a File with no ``md5sum``, the run errors (same rules as ``--reupload``).
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

# Same hosts as igvf_utils.connection / IGVF_MODES (search index matches submission target).
_DEFAULT_PORTAL_API_BY_MODE = {
    "prod": "https://api.data.igvf.org",
    "staging": "https://api.staging.igvf.org",
}

_MATRIX_PRINCIPAL_DIMENSION = "CRISPR guide capture"
_MATRIX_SECONDARY_DIMENSIONS = ["gene expression"]


def _resolve_portal_api_for_mode(
    igvf_mode: Optional[str], portal_api_override: Optional[str]
) -> str:
    """
    API origin for urllib ReferenceFile search; aligns with where iu_register posts.
    Only ``prod`` and ``staging`` are accepted unless ``--portal-api`` overrides.
    """
    if portal_api_override:
        return portal_api_override.rstrip("/")
    mode = (igvf_mode or os.environ.get("IGVF_MODE") or "staging").strip()
    if mode not in _DEFAULT_PORTAL_API_BY_MODE:
        raise ValueError(
            f"Invalid portal mode {mode!r}; use prod or staging (or pass --portal-api URL)."
        )
    try:
        from igvf_utils import IGVF_MODES

        if mode in IGVF_MODES:
            return IGVF_MODES[mode]["url"].rstrip("/")
    except ImportError:
        pass
    return _DEFAULT_PORTAL_API_BY_MODE[mode]


def _submission_igvf_mode(cli_mode: Optional[str]) -> str:
    """Mode string for ``iu_register -m`` (prod or staging)."""
    return (cli_mode or os.environ.get("IGVF_MODE") or "staging").strip()


def _run_iu_register(
    igvf_mode: str,
    profile_id: str,
    tsv_path: str,
    *,
    dry_run: bool,
    no_upload_file: bool,
) -> int:
    """
    Invoke igvf_utils ``iu_register`` (console entry or ``python -m``).
    Returns process exit code.
    """
    cmd: List[str]
    # Prefer `iu_register` on PATH when installed (pip install igvf-utils).
    exe = shutil.which("iu_register")
    if exe:
        cmd = [exe, "-m", igvf_mode, "-p", profile_id, "-i", tsv_path]
    else:
        cmd = [
            sys.executable,
            "-m",
            "igvf_utils.MetaDataRegistration.iu_register",
            "-m",
            igvf_mode,
            "-p",
            profile_id,
            "-i",
            tsv_path,
        ]
    if dry_run:
        cmd.append("--dry-run")
    effective_no_upload = no_upload_file or dry_run
    if effective_no_upload:
        cmd.append("--no-upload-file")
    print("Running:", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=False).returncode


def _run_gsutil_ls(prefix: str) -> List[str]:
    if prefix.startswith("gs://"):
        prefix = prefix.rstrip("/") + "/"
        proc = subprocess.run(
            ["gsutil", "ls", prefix],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"gsutil ls failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}"
            )
        return [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]

    root = Path(prefix).expanduser().resolve()
    if not root.exists():
        raise RuntimeError(f"Input directory does not exist: {str(root)!r}")
    if not root.is_dir():
        raise RuntimeError(f"Input path is not a directory: {str(root)!r}")
    return [str(p) for p in sorted(root.iterdir(), key=lambda x: x.name)]


def _first_n_files(uris: Sequence[str], n: int) -> List[str]:
    files: List[str] = []
    for u in uris:
        if u.startswith("gs://"):
            if not u.endswith("/"):
                files.append(u)
        else:
            p = Path(u)
            if p.exists() and p.is_file():
                files.append(str(p))
    return files[:n]


def _tabular_submission_rank(row: List[str]) -> int:
    """
    Ensure tabular rows are POSTed in derived_from-safe order.
    per_guide records must exist before per_element records.
    """
    alias = row[3].lower() if len(row) > 3 else ""
    if "cis_per_guide" in alias:
        return 0
    if "trans_per_guide" in alias:
        return 1
    if "cis_per_element" in alias:
        return 2
    if "trans_per_element" in alias:
        return 3
    return 99


def _classify(basename: str) -> Tuple[str, str]:
    """
    Returns (profile_id, content_type) for igvf_utils Connection / iu_register.
    """
    lower = basename.lower()
    if lower.endswith(".h5mu") or "inference_mudata" in lower:
        return (
            "matrix_file",
            "annotated multimodal CRISPR matrix",
        )
    if "per_guide" in lower:
        return "tabular_file", "differential guide quantifications"
    if "per_element" in lower:
        return "tabular_file", "differential element quantifications"
    raise ValueError(f"Cannot classify pipeline artifact: {basename!r}")


def _alias_stem_from_basename(basename: str) -> str:
    """
    Return basename without common file extensions used by this pipeline.
    Handles compound extensions like .tsv.gz.
    """
    lower = basename.lower()
    for ext in (".tsv.gz", ".fastq.gz", ".yaml.gz", ".h5mu", ".tsv", ".gz", ".yaml"):
        if lower.endswith(ext):
            return basename[: -len(ext)]
    if "." in basename:
        return basename.rsplit(".", 1)[0]
    return basename


def _md5_and_size_gs(uri: str) -> Tuple[str, int]:
    try:
        from igvf_utils.utils import calculate_file_size, calculate_md5sum
    except ImportError as exc:
        raise RuntimeError(
            "Computing md5sum/size for gs:// requires igvf_utils "
            "(pip install path/to/igvf_utils or set PYTHONPATH)."
        ) from exc
    return calculate_md5sum(uri), calculate_file_size(uri)


def _md5_and_size_local(path: str) -> Tuple[str, int]:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise RuntimeError(f"Local file not found: {str(p)!r}")
    if not p.is_file():
        raise RuntimeError(f"Local path is not a file: {str(p)!r}")

    md5 = hashlib.md5()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            md5.update(chunk)
    return md5.hexdigest(), p.stat().st_size


def _md5_gsutil_hash(uri: str) -> str:
    proc = subprocess.run(
        ["gsutil", "hash", "-m", uri],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"gsutil hash failed: {proc.stderr.strip()}")
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.lower().startswith("hash (md5):"):
            b64 = line.split(":", 1)[1].strip().replace("\t", "")
            # GCS / gsutil report MD5 as base64; portal and igvf_utils expect hex.
            return base64.b64decode(b64).hex()
    raise RuntimeError(f"Could not parse md5 from gsutil hash output:\n{proc.stdout}")


def _md5_for_uri(uri: str, use_gsutil_hash: bool) -> str:
    if not uri.startswith("gs://"):
        md5, _ = _md5_and_size_local(uri)
        return md5
    if use_gsutil_hash:
        return _md5_gsutil_hash(uri)
    try:
        md5, _ = _md5_and_size_gs(uri)
        return md5
    except Exception as exc:
        # igvf_utils GCS lookup can fail when local ADC credentials are expired.
        # Fall back to gsutil, which may still have valid auth context.
        print(
            f"warning: igvf_utils md5 lookup failed for {uri!r}: {exc}; "
            "falling back to gsutil hash -m",
            file=sys.stderr,
        )
        return _md5_gsutil_hash(uri)


def _run_gsutil_cat(uri: str) -> str:
    if not uri.startswith("gs://"):
        p = Path(uri).expanduser().resolve()
        if not p.exists() or not p.is_file():
            raise RuntimeError(f"Local params file not found: {str(p)!r}")
        return p.read_text(encoding="utf-8")
    proc = subprocess.run(
        ["gsutil", "cat", uri],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"gsutil cat failed for {uri!r} ({proc.returncode}): "
            f"{proc.stderr.strip() or proc.stdout.strip()}"
        )
    return proc.stdout


def _search_reference_files_by_md5(api_base: str, md5_hex: str) -> List[str]:
    """
    Query IGVF search for ReferenceFile items with the given md5sum.
    Returns accession strings (e.g. IGVFFI0653VCGH).
    """
    q = urllib.parse.urlencode(
        [("type", "ReferenceFile"), ("md5sum", md5_hex), ("limit", "all")]
    )
    url = urllib.parse.urljoin(api_base.rstrip("/") + "/", "search/?" + q)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.load(resp)
    except urllib.error.HTTPError as exc:
        # IGVF search returns 404 when the query matches no items.
        if exc.code == 404:
            return []
        raise RuntimeError(f"IGVF search failed ({exc.code}): {exc.reason}") from exc
    accessions: List[str] = []
    for item in payload.get("@graph", []):
        acc = item.get("accession")
        if acc and acc not in accessions:
            accessions.append(acc)
    return accessions


def _get_reference_file_object(api_base: str, ref_accession: str) -> Dict[str, Any]:
    """
    Fetch a ReferenceFile object payload from the portal by accession.
    """
    ref = ref_accession.strip()
    url = urllib.parse.urljoin(
        api_base.rstrip("/") + "/",
        f"reference-files/{urllib.parse.quote(ref)}/@@object?format=json",
    )
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.load(resp)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"Failed to fetch ReferenceFile {ref!r} ({exc.code}): {exc.reason}"
        ) from exc
    return payload


def _derived_from_reference(api_base: str, ref_identifier: str) -> Tuple[str, List[str]]:
    ref_obj = _get_reference_file_object(api_base, ref_identifier)
    ref_acc = ref_obj.get("accession") or ref_identifier
    derived_from_raw = ref_obj.get("derived_from") or []
    derived_from_accs: List[str] = []
    for item in derived_from_raw:
        acc = _normalize_link_identifier(item)
        if acc and acc not in derived_from_accs:
            derived_from_accs.append(acc)
    if len(derived_from_accs) < 2:
        raise RuntimeError(
            f"ReferenceFile {ref_acc!r} has insufficient derived_from entries "
            f"for reference_files: {derived_from_raw!r}"
        )
    # Include the index ReferenceFile itself plus upstream references.
    return ref_acc, [ref_acc, *derived_from_accs]


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes", "y"}:
            return True
        if lowered in {"false", "0", "no", "n"}:
            return False
    raise RuntimeError(f"Expected boolean for use_igvf_reference, got {value!r}")


def _run_token_from_params_uri(params_uri: str) -> str:
    """
    Build a deterministic run token from params filename:
    params_YYYY-MM-DD_HH-MM-SS.json -> YYYYMMDD_HHMMSS
    Fallback: short hash of params URI/path.
    """
    name = Path(params_uri).name
    m = re.match(r"^params_(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})\.json$", name)
    if m:
        y, mo, d, h, mi, s = m.groups()
        return f"{y}{mo}{d}_{h}{mi}{s}"
    return hashlib.sha1(params_uri.encode("utf-8")).hexdigest()[:10]


def _detect_use_igvf_reference(gs_prefix: str) -> Tuple[bool, str, str]:
    base = gs_prefix.rstrip("/")
    params_prefix = f"{base}/pipeline_info/" if base.startswith("gs://") else str(Path(base) / "pipeline_info")
    params_candidates = [u for u in _run_gsutil_ls(params_prefix) if not u.endswith("/")]
    params_candidates = [u for u in params_candidates if "params" in Path(u).name and u.endswith(".json")]
    if not params_candidates:
        raise RuntimeError(
            f"Could not find pipeline params JSON under {params_prefix!r} "
            "(expected pipeline_info/params*.json)."
        )
    params_uri = sorted(params_candidates)[-1]
    try:
        payload = json.loads(_run_gsutil_cat(params_uri))
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON in {params_uri!r}: {exc}") from exc
    if "use_igvf_reference" not in payload:
        raise RuntimeError(
            f"'use_igvf_reference' missing in pipeline params file {params_uri!r}"
        )
    run_token = _run_token_from_params_uri(params_uri)
    return _parse_bool(payload["use_igvf_reference"]), params_uri, run_token


def _resolve_reference_files(
    gs_prefix: str,
    api_base: str,
    reference_files_arg: str,
) -> Tuple[List[str], List[str], str]:
    """
    Returns reference_files and log lines.
    Uses pipeline params to determine whether to use the IGVF default reference
    or require explicit --reference-files.
    """
    log: List[str] = []
    use_igvf_reference, params_uri, run_token = _detect_use_igvf_reference(gs_prefix)
    log.append(f"pipeline params: {params_uri}")
    log.append(f"use_igvf_reference={use_igvf_reference}")
    log.append(f"run_token={run_token}")

    if use_igvf_reference:
        ref_seed = reference_files_arg.strip() if reference_files_arg.strip() else "IGVFFI9561BASO"
        ref_acc, derived_from_accs = _derived_from_reference(api_base, ref_seed)
        log.append(
            f"use_igvf_reference=true -> ReferenceFile {ref_acc} derived_from -> "
            f"{','.join(derived_from_accs)}"
        )
        return derived_from_accs, log, run_token

    raw = [x.strip() for x in reference_files_arg.split(",") if x.strip()]
    if not raw:
        raise RuntimeError(
            "Pipeline params set use_igvf_reference=false; --reference-files is required."
        )
    normalized = [_normalize_link_identifier(x) for x in raw]
    out: List[str] = []
    for item in normalized:
        if item and item not in out:
            out.append(item)
    if not out:
        raise RuntimeError(
            "No valid values parsed from --reference-files; provide comma-separated accessions/aliases."
        )
    log.append(f"use_igvf_reference=false -> using --reference-files: {','.join(out)}")
    return out, log, run_token


# --- Re-upload (same analysis set, md5-matched File records) -----------------


def _file_set_accession_from_embed(file_rec: dict) -> Optional[str]:
    fs = file_rec.get("file_set")
    if isinstance(fs, dict):
        return fs.get("accession")
    return None


def _ensure_file_set_accession(conn: Any, file_rec: dict) -> str:
    acc = _file_set_accession_from_embed(file_rec)
    if acc:
        return acc
    fid = file_rec.get("accession")
    if not fid:
        href = file_rec.get("@id", "")
        parts = [p for p in href.strip("/").split("/") if p]
        if len(parts) >= 2:
            fid = parts[-1]
    if not fid:
        return ""
    full = conn.get(rec_ids=fid, ignore404=False, database=True, frame="object")
    return _file_set_accession_from_embed(full) or ""


def _resolve_file_set_accession(conn: Any, ref: str) -> str:
    rec = conn.get(rec_ids=ref.strip(), ignore404=False, database=True, frame="object")
    acc = rec.get("accession")
    if not acc:
        raise SystemExit(f"Could not resolve accession for file set {ref!r}")
    return acc


def _normalize_link_identifier(value: Any) -> str:
    """
    Normalize portal link-like values to identifier tails.
    Examples:
      - "/awards/HG012053/" -> "HG012053"
      - "/labs/charles-gersbach/" -> "charles-gersbach"
      - "HG012053" -> "HG012053"
    """
    if isinstance(value, dict):
        if value.get("accession"):
            return str(value["accession"]).strip()
        value = value.get("@id") or value.get("uuid") or ""
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("/"):
        parts = [p for p in text.strip("/").split("/") if p]
        if parts:
            return parts[-1]
    return text


def _resolve_award_and_lab(conn: Any, file_set_ref: str) -> Tuple[str, str]:
    """
    Return (award, lab) from analysis-set metadata.
    """
    rec = conn.get(rec_ids=file_set_ref.strip(), ignore404=False, database=True, frame="object")
    award = _normalize_link_identifier(rec.get("award"))
    lab = _normalize_link_identifier(rec.get("lab"))

    missing = []
    if not award:
        missing.append("award")
    if not lab:
        missing.append("lab")
    if missing:
        raise SystemExit(
            f"Could not resolve {', '.join(missing)} from analysis set {file_set_ref!r}."
        )
    return award, lab


def _patch_uniform_pipeline_status_completed(conn: Any, file_set_ref: str, *, dry_run: bool) -> None:
    """
    PATCH the target analysis set with uniform_pipeline_status=completed.
    """
    rec = conn.get(rec_ids=file_set_ref.strip(), ignore404=False, database=True, frame="object")
    rec_id = rec.get("@id")
    if not rec_id:
        raise SystemExit(
            f"Could not determine @id for analysis set {file_set_ref!r} for status patch."
        )
    if dry_run:
        print(
            f"[dry-run] Would patch {rec_id} with uniform_pipeline_status='completed'.",
            flush=True,
        )
        return
    conn.set_submission(True)
    conn.patch(
        {
            conn.IGVFID_KEY: rec_id,
            "uniform_pipeline_status": "completed",
        },
        extend_array_values=False,
    )
    print(
        f"Patched {rec_id} with uniform_pipeline_status='completed'.",
        flush=True,
    )


def _search_files_by_md5(conn: Any, md5_hex: str) -> List[dict]:
    return conn.search([("type", "File"), ("md5sum", md5_hex)])


def _list_files_in_file_set(conn: Any, file_set_accession: str) -> List[dict]:
    return conn.search([("type", "File"), ("file_set.accession", file_set_accession)])


def _index_md5_for_file_set(
    conn: Any, file_set_accession: str
) -> Tuple[Dict[str, str], List[str], List[str]]:
    """
    Returns (md5_hex -> file accession, log lines, file accessions missing md5sum).
    """
    rows = _list_files_in_file_set(conn, file_set_accession)
    md5_map: Dict[str, str] = {}
    log_lines: List[str] = []
    missing_md5_accessions: List[str] = []
    for rec in rows:
        acc = rec.get("accession")
        md5 = rec.get("md5sum")
        if not md5:
            full = conn.get(rec_ids=acc, ignore404=False, database=True, frame="object")
            md5 = full.get("md5sum")
        if not md5:
            if acc:
                missing_md5_accessions.append(acc)
            log_lines.append(f"  skip indexing {acc}: no md5sum on portal")
            continue
        if md5 in md5_map and md5_map[md5] != acc:
            raise SystemExit(
                f"Ambiguous metadata on file set {file_set_accession!r}: md5sum {md5} "
                f"on both {md5_map[md5]!r} and {acc!r}"
            )
        md5_map[md5] = acc
    return md5_map, log_lines, missing_md5_accessions


def _partition_artifacts_for_post_and_upload(
    conn: Any,
    on_set_index: Dict[str, str],
    artifacts: List[Tuple[str, str, str, List[str]]],
) -> Tuple[List[List[str]], List[List[str]], List[Tuple[str, str]], List[str]]:
    """
    For each artifact (uri, md5, profile, row): if md5 is not on the file set, row goes to POST;
    if md5 is on the set and upload_status is not a terminal state (``validated`` or
    ``validation exempted``), queue upload; otherwise skip (log only).

    Returns (post_matrix_rows, post_tabular_rows, upload_jobs [(accession, uri)], skip_msgs).
    """
    post_matrix: List[List[str]] = []
    post_tabular: List[List[str]] = []
    uploads: List[Tuple[str, str]] = []
    skip_msgs: List[str] = []

    for uri, md5, profile, row in artifacts:
        if md5 not in on_set_index:
            if profile == "matrix_file":
                post_matrix.append(row)
            else:
                post_tabular.append(row)
            continue

        acc = on_set_index[md5]
        rec = conn.get(rec_ids=acc, ignore404=False, database=True, frame="object")
        st = rec.get("upload_status")
        if st in ["validated", "validation exempted"]:
            skip_msgs.append(
                f"skip POST/upload: {acc} upload_status={st!r} (no re-upload; md5={md5})"
            )
        else:
            uploads.append((acc, uri))
            skip_msgs.append(
                f"will upload blob for {acc} (upload_status={st!r}, md5={md5})"
            )

    return post_matrix, post_tabular, uploads, skip_msgs


def _validate_file_set_md5s_match_sources(
    on_set_index: Dict[str, str],
    source_md5_set: Set[str],
    target_fs_acc: str,
) -> None:
    """
    Every File on the analysis set must have an md5sum that appears among the
    pipeline output blobs being re-uploaded.
    """
    orphan = []
    for md5_hex, file_acc in on_set_index.items():
        if md5_hex not in source_md5_set:
            orphan.append(f"{file_acc} (md5={md5_hex})")
    if orphan:
        raise SystemExit(
            f"File(s) on analysis set {target_fs_acc!r} do not match any pipeline output md5: "
            + "; ".join(orphan)
        )


def _pick_reupload_targets_for_source(
    conn: Any,
    source_md5: str,
    target_fs_acc: str,
    on_set_index: Dict[str, str],
) -> List[str]:
    hits = _search_files_by_md5(conn, source_md5)
    on_target: List[str] = []
    on_other: List[Tuple[str, str]] = []

    for rec in hits:
        facc = rec.get("accession")
        fs_acc = _file_set_accession_from_embed(rec) or _ensure_file_set_accession(conn, rec)
        if not fs_acc:
            raise SystemExit(
                f"Could not determine file_set for file {facc!r} (md5={source_md5})"
            )
        if fs_acc == target_fs_acc:
            if facc and facc not in on_target:
                on_target.append(facc)
        else:
            on_other.append((facc, fs_acc))

    if on_other:
        detail = ", ".join(f"{a} (file_set {fs})" for a, fs in on_other)
        raise SystemExit(
            f"md5sum {source_md5} matches file(s) outside file set {target_fs_acc!r}: {detail}. "
            "Resolve before re-uploading."
        )

    if len(on_target) > 1:
        raise SystemExit(
            f"md5sum {source_md5} matches multiple files on file set {target_fs_acc!r}: {on_target!r}"
        )

    idx_acc = on_set_index.get(source_md5)
    if len(on_target) == 1:
        acc = on_target[0]
        if idx_acc and idx_acc != acc:
            raise SystemExit(
                f"Inconsistent portal data for md5sum {source_md5}: file_set listing has {idx_acc!r}, "
                f"md5 search has {acc!r}"
            )
        return [acc]

    if idx_acc:
        return [idx_acc]

    return []


def _md5_source_uri(uri: str, use_gsutil_hash: bool) -> str:
    return _md5_for_uri(uri, use_gsutil_hash)


def _main_reupload(argv: Sequence[str]) -> int:
    try:
        from igvf_utils.parent_argparser import igvf_login_parser
    except ImportError as exc:
        raise SystemExit(
            "The --reupload mode requires igvf_utils (e.g. pip install / path to igvf_utils). "
            f"Original error: {exc}"
        ) from exc

    parser = argparse.ArgumentParser(
        description="Re-upload GCS (or local) objects onto existing File records on an analysis set.",
        parents=[igvf_login_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python accession_outputs.py --reupload -m staging --analysis-set IGVFDS5057HJKP "
            "gs://my-bucket/pipeline_runs/run_01/\n"
            "  python accession_outputs.py --reupload -m staging --analysis-set IGVFDS5057HJKP "
            "--source-uri gs://my-bucket/pipeline_runs/run_01/inference_mudata.h5mu "
            "--source-uri gs://my-bucket/pipeline_runs/run_01/per_guide_results.tsv.gz"
        ),
    )
    parser.add_argument(
        "--analysis-set",
        required=True,
        dest="analysis_set",
        help="Target analysis set accession or alias (e.g. IGVFDS5057HJKP).",
    )
    parser.add_argument(
        "gs_prefix",
        nargs="?",
        default=None,
        help="Run directory gs:// prefix (optional if --source-uri is used).",
    )
    parser.add_argument(
        "--source-uri",
        action="append",
        dest="source_uris",
        default=None,
        help="Explicit gs:// or local path; repeatable. Overrides default five files from prefix.",
    )
    parser.add_argument(
        "--use-gsutil-hash",
        action="store_true",
        help="Compute GCS md5 via gsutil hash -m (default: igvf_utils).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not upload; only print planned actions.",
    )
    parser.add_argument(
        "--no-log-file",
        action="store_true",
        help="Disable igvf_utils Connection file logging.",
    )
    args = parser.parse_args(list(argv))

    try:
        from igvf_utils.connection import Connection
    except ImportError as exc:
        raise SystemExit(
            "Could not import igvf_utils.connection.Connection for --reupload. "
            f"Original error: {exc}"
        ) from exc

    if args.source_uris:
        sources = list(args.source_uris)
    else:
        if not args.gs_prefix:
            parser.error("gs_prefix is required unless --source-uri is given")
        all_uris = _run_gsutil_ls(args.gs_prefix)
        sources = _first_n_files(all_uris, 5)
        if len(sources) < 5:
            parser.error(
                f"Expected at least 5 objects under {args.gs_prefix!r}; found {len(sources)}"
            )

    conn = Connection(
        igvf_mode=args.igvf_mode,
        dry_run=args.dry_run,
        submission=False,
        no_log_file=args.no_log_file,
    )
    target_fs = _resolve_file_set_accession(conn, args.analysis_set)

    on_set_index, index_logs, missing_md5_accs = _index_md5_for_file_set(conn, target_fs)
    print(f"Target file_set accession: {target_fs}")
    print(f"Indexed {len(on_set_index)} md5 value(s) from Files on this set.")
    for line in index_logs:
        print(line)

    source_md5_set: Set[str] = set()
    source_pairs: List[Tuple[str, str]] = []
    for uri in sources:
        try:
            md5 = _md5_source_uri(uri, args.use_gsutil_hash)
        except Exception as exc:
            print(f"error: could not hash {uri!r}: {exc}", file=sys.stderr)
            return 1
        source_md5_set.add(md5)
        source_pairs.append((uri, md5))

    if missing_md5_accs:
        print(
            f"error: cannot validate File(s) with no md5sum on the portal: {missing_md5_accs!r}",
            file=sys.stderr,
        )
        return 1

    try:
        _validate_file_set_md5s_match_sources(on_set_index, source_md5_set, target_fs)
    except SystemExit as exc:
        print(f"error: {exc.args[0]}", file=sys.stderr)
        return 1

    conn.set_submission(True)

    exit_code = 0
    for uri, md5 in source_pairs:
        try:
            targets = _pick_reupload_targets_for_source(conn, md5, target_fs, on_set_index)
        except SystemExit as exc:
            print(f"error: {exc.args[0]}", file=sys.stderr)
            exit_code = 1
            continue

        if not targets:
            print(
                f"skip: no File with md5sum={md5} on file set {target_fs!r} (source {uri!r})"
            )
            continue

        file_id = targets[0]
        print(f"reupload {uri!r} -> {file_id} (md5={md5})")
        try:
            conn.upload_file(file_id, file_path=uri, set_md5sum=True)
        except Exception as exc:
            print(f"error: upload failed for {file_id}: {exc}", file=sys.stderr)
            exit_code = 1

    return exit_code


def _main_generate(argv: Sequence[str]) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python accession_outputs.py -m staging --analysis-set IGVFDS5057HJKP "
            "gs://my-bucket/pipeline_runs/run_01/\n"
            "  python accession_outputs.py -m staging --analysis-set charles-gersbach:my_analysis_set "
            "gs://my-bucket/pipeline_runs/run_01/\n\n"
            "  python accession_outputs.py -m staging --analysis-set IGVFDS5057HJKP "
            "--reference-files IGVFFI1234ABCD,IGVFFI5678EFGH "
            "gs://my-bucket/pipeline_runs/run_01/\n\n"
            "  python accession_outputs.py -m staging --analysis-set IGVFDS5057HJKP "
            "/path/to/local/run_01/\n\n"
            "Tip: add --reupload to push blobs onto existing File records only."
        ),
    )
    p.add_argument(
        "-m",
        "--igvf-mode",
        default=None,
        dest="igvf_mode",
        help="Portal for ReferenceFile lookup and for ``iu_register -m`` (prod or staging only). "
        "Default: IGVF_MODE env var, else staging.",
    )
    p.add_argument(
        "gs_prefix",
        help="Input run directory prefix: gs://bucket/.../Engreitz_ccPerturb or /local/path/run",
    )
    p.add_argument(
        "--analysis-set",
        required=True,
        dest="analysis_set",
        help="analysis set alias or accession (e.g. charles-gersbach:my_analysis_set or IGVFDS5057HJKP).",
    )
    p.add_argument(
        "--portal-api",
        default=None,
        dest="portal_api",
        help="Override API origin for ReferenceFile search (default: from -m / IGVF_MODE, "
        "same host as igvf_utils submissions).",
    )
    p.add_argument(
        "--reference-files",
        default="",
        help="Optional comma-separated reference file accessions/aliases. "
        "If pipeline params set use_igvf_reference=true, defaults to IGVFFI9561BASO "
        "(resolved to itself plus its derived_from values). If use_igvf_reference=false, this argument is required.",
    )
    p.add_argument(
        "--use-gsutil-hash",
        action="store_true",
        help="Use gsutil hash -m for md5 instead of igvf_utils.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="No portal writes: skip uploads and pass --dry-run to iu_register.",
    )
    p.add_argument(
        "--no-upload-file",
        action="store_true",
        help="POST metadata only: skip automatic blob uploads and pass --no-upload-file to iu_register.",
    )
    args = p.parse_args(list(argv) if argv is not None else None)

    try:
        portal_api = _resolve_portal_api_for_mode(args.igvf_mode, args.portal_api)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    try:
        ref_files, ref_log, run_token = _resolve_reference_files(
            args.gs_prefix,
            portal_api,
            args.reference_files,
        )
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    all_uris = _run_gsutil_ls(args.gs_prefix)
    picked = _first_n_files(all_uris, 5)
    if len(picked) < 5:
        raise SystemExit(
            f"Expected at least 5 files under {args.gs_prefix!r}; found {len(picked)}: {picked}"
        )

    submit_mode = _submission_igvf_mode(args.igvf_mode)
    if submit_mode not in _DEFAULT_PORTAL_API_BY_MODE:
        print(
            f"error: submission requires igvf mode prod or staging (got {submit_mode!r}); "
            "set -m prod or -m staging (or IGVF_MODE).",
            file=sys.stderr,
        )
        return 1

    print("\nreference_files resolution:")
    print(f"  ReferenceFile md5 search API: {portal_api}")
    for line in ref_log:
        print(f"  {line}")
    print(f"  matrix_file & tabular_file reference_files -> {','.join(ref_files)}")

    try:
        from igvf_utils.connection import Connection
    except ImportError as exc:
        print(
            f"error: portal sync requires igvf_utils (Connection). {exc}",
            file=sys.stderr,
        )
        return 1

    conn = Connection(
        igvf_mode=submit_mode,
        dry_run=args.dry_run,
        submission=False,
        no_log_file=False,
    )
    target_fs = _resolve_file_set_accession(conn, args.analysis_set)
    try:
        award, lab = _resolve_award_and_lab(conn, args.analysis_set)
    except SystemExit as exc:
        print(f"error: {exc.args[0]}", file=sys.stderr)
        return 1
    print(f"Using award={award!r} and lab={lab!r} for submitted rows.")
    alias_prefix = f"igvf:{target_fs}_uniform_perturb_seq_pipeline_{run_token}_"
    print(f"Using alias prefix: {alias_prefix!r}")

    # Build artifact metadata first so intra-run derived_from links can be resolved.
    built: List[Dict[str, str]] = []

    for uri in picked:
        base = uri.rsplit("/", 1)[-1]
        profile, ctype = _classify(base)
        md5 = _md5_for_uri(uri, args.use_gsutil_hash)

        alias_stem = _alias_stem_from_basename(base)
        alias = f"{alias_prefix}{re.sub(r'[^0-9a-zA-Z._-]+', '_', alias_stem)}"

        built.append(
            {
                "uri": uri,
                "base": base,
                "profile": profile,
                "content_type": ctype,
                "md5": md5,
                "alias": alias,
            }
        )

    matrix_aliases = [x["alias"] for x in built if x["profile"] == "matrix_file"]
    matrix_alias = matrix_aliases[0] if matrix_aliases else ""
    cis_guide_alias = next(
        (
            x["alias"]
            for x in built
            if x["profile"] == "tabular_file" and "cis_per_guide" in x["base"].lower()
        ),
        "",
    )
    trans_guide_alias = next(
        (
            x["alias"]
            for x in built
            if x["profile"] == "tabular_file" and "trans_per_guide" in x["base"].lower()
        ),
        "",
    )

    # (uri, md5, profile_id, tsv_row_cells)
    artifacts: List[Tuple[str, str, str, List[str]]] = []
    for item in built:
        uri = item["uri"]
        md5 = item["md5"]
        profile = item["profile"]
        ctype = item["content_type"]
        alias = item["alias"]
        base_lower = item["base"].lower()

        if profile == "matrix_file":
            row = [
                award,
                lab,
                md5,
                alias,
                uri,
                "h5mu",
                ctype,
                args.analysis_set,
                ",".join(ref_files),
                _MATRIX_PRINCIPAL_DIMENSION,
                ",".join(_MATRIX_SECONDARY_DIMENSIONS),
                "",
            ]
        else:
            derived_from = ""
            if "per_guide" in base_lower:
                derived_from = matrix_alias
            elif "cis_per_element" in base_lower:
                derived_from = cis_guide_alias
            elif "trans_per_element" in base_lower:
                derived_from = trans_guide_alias
            row = [
                award,
                lab,
                md5,
                alias,
                uri,
                "tsv",
                ctype,
                args.analysis_set,
                ",".join(ref_files),
                "false",
                derived_from,
            ]
        artifacts.append((uri, md5, profile, row))

    n_matrix = sum(1 for _, _, prof, _ in artifacts if prof == "matrix_file")
    n_tabular = sum(1 for _, _, prof, _ in artifacts if prof == "tabular_file")
    if n_matrix != 1 or n_tabular != 4:
        raise SystemExit(
            f"Expected 1 matrix_file and 4 tabular_file rows; got matrix={n_matrix}, tabular={n_tabular}. "
            "Check the first five objects under the prefix match inference_mudata.h5mu and four *results.tsv.gz files."
        )
    if not matrix_alias:
        raise SystemExit("Could not resolve matrix_file alias for derived_from links.")
    has_cis_element = any("cis_per_element" in x["base"].lower() for x in built)
    has_trans_element = any("trans_per_element" in x["base"].lower() for x in built)
    if has_cis_element and not cis_guide_alias:
        raise SystemExit(
            "Found cis_per_element output but could not resolve cis_per_guide alias for derived_from."
        )
    if has_trans_element and not trans_guide_alias:
        raise SystemExit(
            "Found trans_per_element output but could not resolve trans_per_guide alias for derived_from."
        )

    matrix_header = "\t".join(
        [
            "award",
            "lab",
            "md5sum",
            "aliases",
            "submitted_file_name",
            "file_format",
            "content_type",
            "file_set",
            "reference_files",
            "principal_dimension",
            "secondary_dimensions",
            "derived_from",
        ]
    )
    tab_header = "\t".join(
        [
            "award",
            "lab",
            "md5sum",
            "aliases",
            "submitted_file_name",
            "file_format",
            "content_type",
            "file_set",
            "reference_files",
            "controlled_access",
            "derived_from",
        ]
    )

    on_set_index, index_logs, missing_md5_accs = _index_md5_for_file_set(conn, target_fs)
    source_md5_set: Set[str] = {md5 for _, md5, _, _ in artifacts}

    print(f"\nAnalysis set {target_fs!r}: {len(on_set_index)} indexed File(s) with md5sum.")
    for line in index_logs:
        print(line)

    if missing_md5_accs:
        print(
            f"error: cannot validate File(s) with no md5sum on the portal: {missing_md5_accs!r}",
            file=sys.stderr,
        )
        return 1

    if on_set_index:
        try:
            _validate_file_set_md5s_match_sources(on_set_index, source_md5_set, target_fs)
        except SystemExit as exc:
            print(f"error: {exc.args[0]}", file=sys.stderr)
            return 1

    post_matrix_rows, post_tabular_rows, upload_jobs, sync_msgs = (
        _partition_artifacts_for_post_and_upload(conn, on_set_index, artifacts)
    )
    post_tabular_rows.sort(key=_tabular_submission_rank)

    print("\nPortal sync (existing md5 on analysis set):")
    for msg in sync_msgs:
        print(f"  {msg}")

    if upload_jobs and not args.no_upload_file:
        conn.set_submission(True)
        print(
            f"\nUploading {len(upload_jobs)} file(s) "
            f"(upload_status neither validated nor validation exempted)...",
            flush=True,
        )
        for acc, uri in upload_jobs:
            print(f"  upload {uri!r} -> {acc}", flush=True)
            try:
                conn.upload_file(acc, file_path=uri, set_md5sum=True)
            except Exception as exc:
                print(f"error: upload failed for {acc}: {exc}", file=sys.stderr)
                return 1
    elif upload_jobs and args.no_upload_file:
        print(
            "\nSkipping blob uploads (--no-upload-file); existing records still need validation.",
            flush=True,
        )

    if not post_matrix_rows and not post_tabular_rows:
        print("\nNo new File records to POST (all md5s already on analysis set).", flush=True)
        try:
            _patch_uniform_pipeline_status_completed(conn, args.analysis_set, dry_run=args.dry_run)
        except Exception as exc:
            print(
                f"error: failed to patch uniform_pipeline_status on {args.analysis_set!r}: {exc}",
                file=sys.stderr,
            )
            return 1
        return 0

    with tempfile.TemporaryDirectory(prefix="igvf_register_") as tmp:
        rc = 0
        if post_matrix_rows:
            mpath = os.path.join(tmp, "register_matrix_file.tsv")
            with open(mpath, "w", encoding="utf-8") as fh:
                fh.write(matrix_header + "\n")
                for r in post_matrix_rows:
                    fh.write("\t".join(r) + "\n")
            print(
                f"\nSubmitting {len(post_matrix_rows)} matrix_file row(s) via iu_register...",
                flush=True,
            )
            rc = _run_iu_register(
                submit_mode,
                "matrix_file",
                mpath,
                dry_run=args.dry_run,
                no_upload_file=args.no_upload_file,
            )
            if rc != 0:
                return rc

        if post_tabular_rows:
            tpath = os.path.join(tmp, "register_tabular_file.tsv")
            with open(tpath, "w", encoding="utf-8") as fh:
                fh.write(tab_header + "\n")
                for r in post_tabular_rows:
                    fh.write("\t".join(r) + "\n")
            print(
                f"\nSubmitting {len(post_tabular_rows)} tabular_file row(s) via iu_register...",
                flush=True,
            )
            rc = _run_iu_register(
                submit_mode,
                "tabular_file",
                tpath,
                dry_run=args.dry_run,
                no_upload_file=args.no_upload_file,
            )
        if rc != 0:
            return rc
        try:
            _patch_uniform_pipeline_status_completed(conn, args.analysis_set, dry_run=args.dry_run)
        except Exception as exc:
            print(
                f"error: failed to patch uniform_pipeline_status on {args.analysis_set!r}: {exc}",
                file=sys.stderr,
            )
            return 1
        return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--reupload", action="store_true")
    pre_args, rest = pre.parse_known_args(argv_list)
    if pre_args.reupload:
        return _main_reupload(rest)
    return _main_generate(rest)


if __name__ == "__main__":
    raise SystemExit(main())
