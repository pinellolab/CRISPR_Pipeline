#!/usr/bin/env python3
"""
Generate, submit, and monitor six Nextflow QC sweep runs.

Default command shape per run:
  NXF_VER=26.03.2-edge nextflow run main.nf -profile google \
    -c <generated_config> \
    --input ../sample_metadata_gcp_2026_02_26_patched.csv \
    --outdir gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/cc_olga_version_runs/jesse_{info} \
    --with_tower

The generated configs include the base Nextflow config, then override the run
specific params. Run this script from anywhere; Nextflow itself is launched from
the CRISPR_Pipeline directory. The script stores local stdout logs, an event
log, a state file, and a manifest so launched runs can be found again even if
the TUI terminal closes. Fresh launches do not use Nextflow resume; interrupted
runs can be resumed from the TUI with r. Completed runs can be force-relaunched
fresh with f or --force.
"""

from __future__ import annotations

import argparse
import curses
import datetime as dt
import json
import os
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent

DEFAULT_NXF_VER = "26.03.2-edge"
DEFAULT_BASE_CONFIG = "nextflow_cc.config"
DEFAULT_CONFIG_DIR = "run_configs/qc_sweep_jesse"
DEFAULT_LOG_DIR = "run_logs/qc_sweep_jesse"
DEFAULT_STATE_FILE = ".qc_sweep_tui_state.json"
DEFAULT_EVENT_LOG = "run_logs/qc_sweep_jesse/qc_sweep_events.jsonl"
DEFAULT_MANIFEST_FILE = "run_logs/qc_sweep_jesse/runs_manifest.tsv"
DEFAULT_INPUT = "../sample_metadata_gcp_2026_02_26_patched.csv"
DEFAULT_OUTDIR_TEMPLATE = (
    "gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/"
    "cc_olga_version_runs/jesse_{info}"
)
DEFAULT_RESUME_ID = None
DEFAULT_WITH_TOWER = True
DEFAULT_TOWER_ACCESS_TOKEN_ENV = "TOWER_ACCESS_TOKEN"

RUNNING = "running"
PENDING = "pending"
SUCCEEDED = "succeeded"
INTERRUPTED = "interrupted"
DRY_RUN = "dry-run"
RESUMABLE_STATUSES = {INTERRUPTED}
NEXTFLOW_SESSION_RE = re.compile(
    r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b"
)
NEXTFLOW_STATUS_TO_STATE = {
    "OK": SUCCEEDED,
    "ERR": INTERRUPTED,
    "-": RUNNING,
}


@dataclass(frozen=True)
class RunSpec:
    info: str
    qc_barcode_filter: str
    qc_min_genes_per_cell: Optional[int]

    @property
    def min_genes_display(self) -> str:
        if self.qc_min_genes_per_cell is None:
            return "ignored"
        return str(self.qc_min_genes_per_cell)


RUN_SPECS: Tuple[RunSpec, ...] = (
    RunSpec("min_genes_200_filter_none", "none", 200),
    RunSpec("min_genes_500_filter_none", "none", 500),
    RunSpec("min_genes_800_filter_none", "none", 800),
    RunSpec("min_genes_2000_filter_none", "none", 2000),
    RunSpec("barcode_filter_knee", "knee", None),
    RunSpec("barcode_filter_knee2", "knee2", None),
)


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def resolve_in_pipeline(path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    return SCRIPT_DIR / path


def groovy_quote(value: str) -> str:
    return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"


def read_tail(path: Path, max_bytes: int = 12000) -> str:
    try:
        size = path.stat().st_size
        with path.open("rb") as handle:
            if size > max_bytes:
                handle.seek(size - max_bytes)
            data = handle.read()
    except FileNotFoundError:
        return ""
    except OSError as exc:
        return f"<could not read log: {exc}>"
    return data.decode("utf-8", errors="replace")


def tail_lines(path: Path, limit: int) -> List[str]:
    text = read_tail(path)
    if not text:
        return []
    return text.splitlines()[-limit:]


def infer_status_from_log(path: Path) -> Optional[str]:
    text = read_tail(path, max_bytes=50000).lower()
    if not text:
        return None

    success_markers = (
        "pipeline completed successfully",
        "completed successfully",
        "execution complete",
        "workflow completed",
    )
    if any(marker in text for marker in success_markers):
        return SUCCEEDED

    interrupted_markers = (
        "execution cancelled",
        "keyboardinterrupt",
        "interrupted",
        "session aborted",
        "killing running tasks",
        "terminated",
        "error ~",
    )
    if any(marker in text for marker in interrupted_markers):
        return INTERRUPTED

    return None


def pid_is_alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def discover_local_nextflow_processes(args: argparse.Namespace) -> Dict[str, Dict[str, object]]:
    try:
        completed = subprocess.run(
            ["ps", "-Ao", "pid=,command="],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return {}
    if completed.returncode != 0:
        return {}

    processes: Dict[str, Dict[str, object]] = {}
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        pid_text, _, command = stripped.partition(" ")
        if "nextflow" not in command or " run main.nf " not in f" {command} ":
            continue
        info = history_info_from_command(command, args)
        if info is None:
            continue
        try:
            pid = int(pid_text)
        except ValueError:
            continue
        processes[info] = {"pid": pid, "command": command}
    return processes


def spec_by_info(info: str) -> RunSpec:
    for spec in RUN_SPECS:
        if spec.info == info:
            return spec
    valid = ", ".join(spec.info for spec in RUN_SPECS)
    raise ValueError(f"unknown run info {info!r}; valid values: {valid}")


def config_path_for(spec: RunSpec, config_dir: Path) -> Path:
    return config_dir / f"{spec.info}.config"


def log_path_for(spec: RunSpec, log_dir: Path) -> Path:
    return log_dir / f"{spec.info}.log"


def event_log_path(args: argparse.Namespace) -> Path:
    return resolve_in_pipeline(args.event_log)


def manifest_path(args: argparse.Namespace) -> Path:
    return resolve_in_pipeline(args.manifest_file)


def append_event(
    args: argparse.Namespace,
    event: str,
    spec: Optional[RunSpec] = None,
    record: Optional[Dict[str, object]] = None,
    extra: Optional[Dict[str, object]] = None,
) -> None:
    path = event_log_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, object] = {
        "timestamp": now_iso(),
        "event": event,
        "launcher_pid": os.getpid(),
    }
    if spec is not None:
        payload.update(
            {
                "info": spec.info,
                "qc_barcode_filter": spec.qc_barcode_filter,
                "qc_min_genes_per_cell": spec.qc_min_genes_per_cell,
            }
        )
    if record:
        for key in (
            "status",
            "pid",
            "attempts",
            "returncode",
            "used_resume",
            "resume_id",
            "forced_relaunch",
            "started_at",
            "ended_at",
            "command",
            "config_path",
            "log_path",
            "outdir",
        ):
            if key in record:
                payload[key] = record[key]
    if extra:
        payload.update(extra)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def write_manifest(state: Dict[str, object], args: argparse.Namespace) -> None:
    path = manifest_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    runs = state.get("runs", {})
    if not isinstance(runs, dict):
        runs = {}

    fields = [
        "info",
        "status",
        "pid",
        "attempts",
        "qc_barcode_filter",
        "qc_min_genes_per_cell",
        "started_at",
        "ended_at",
        "returncode",
        "used_resume",
        "resume_id",
        "forced_relaunch",
        "state_path",
        "event_log_path",
        "manifest_path",
        "config_path",
        "log_path",
        "outdir",
        "command",
    ]
    lines = ["\t".join(fields)]
    for spec in RUN_SPECS:
        record = runs.get(spec.info, {})
        if not isinstance(record, dict):
            record = {}
        values = []
        for field in fields:
            if field == "info":
                value = spec.info
            elif field == "qc_barcode_filter":
                value = spec.qc_barcode_filter
            elif field == "qc_min_genes_per_cell":
                value = spec.qc_min_genes_per_cell
            else:
                value = record.get(field, "")
            values.append("" if value is None else str(value).replace("\t", " "))
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_config(spec: RunSpec, base_config: Path) -> str:
    lines = [
        "// Auto-generated by run_qc_sweep_tui.py",
        f"// Generated: {now_iso()}",
        f"// info: {spec.info}",
        f"includeConfig {groovy_quote(str(base_config))}",
        "",
        "params {",
        f"    info = {groovy_quote(spec.info)}",
        f"    QC_barcode_filter = {groovy_quote(spec.qc_barcode_filter)}",
    ]

    if spec.qc_min_genes_per_cell is None:
        lines.append(
            "    // QC_min_genes_per_cell is intentionally not overridden; "
            "knee filters ignore it."
        )
    else:
        lines.append(f"    QC_min_genes_per_cell = {spec.qc_min_genes_per_cell}")

    lines.extend(["}", ""])
    return "\n".join(lines)


def generate_configs(args: argparse.Namespace) -> List[Path]:
    base_config = resolve_in_pipeline(args.base_config)
    config_dir = resolve_in_pipeline(args.config_dir)
    if not base_config.exists():
        raise SystemExit(f"Base config does not exist: {base_config}")

    config_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for spec in RUN_SPECS:
        path = config_path_for(spec, config_dir)
        path.write_text(render_config(spec, base_config), encoding="utf-8")
        written.append(path)
    return written


def outdir_for(spec: RunSpec, args: argparse.Namespace) -> str:
    if "{info}" not in args.outdir_template:
        raise SystemExit("--outdir-template must contain the literal field {info}")
    return args.outdir_template.format(info=spec.info)


def build_command(
    spec: RunSpec,
    args: argparse.Namespace,
    *,
    resume: bool = False,
    resume_id: Optional[str] = None,
) -> List[str]:
    config_dir = resolve_in_pipeline(args.config_dir)
    command = [
        "nextflow",
        "run",
        "main.nf",
        "-profile",
        args.profile,
        "-c",
        str(config_path_for(spec, config_dir)),
        "--input",
        args.input,
        "--outdir",
        outdir_for(spec, args),
    ]
    if args.with_tower:
        command.append("--with_tower")
    if resume:
        command.append("-resume")
        if resume_id:
            command.append(resume_id)
    return command


def display_command(command: List[str], args: argparse.Namespace) -> str:
    return f"NXF_VER={shlex.quote(args.nxf_ver)} {shlex.join(command)}"


def configure_tower_env(env: Dict[str, str], args: argparse.Namespace) -> List[str]:
    """Populate Tower/Seqera env vars without exposing secret values in logs."""
    if not args.with_tower:
        return ["disabled"]

    notes = []
    token_file = getattr(args, "tower_access_token_file", None)
    token_env_name = getattr(args, "tower_access_token_env", DEFAULT_TOWER_ACCESS_TOKEN_ENV)

    if token_file:
        token_path = Path(str(token_file)).expanduser()
        try:
            token = token_path.read_text(encoding="utf-8").strip()
        except OSError as exc:
            raise RuntimeError(f"could not read --tower-access-token-file {token_path}: {exc}") from exc
        if not token:
            raise RuntimeError(f"--tower-access-token-file is empty: {token_path}")
        env["TOWER_ACCESS_TOKEN"] = token
        notes.append("TOWER_ACCESS_TOKEN=file")
    elif token_env_name and os.environ.get(token_env_name):
        env["TOWER_ACCESS_TOKEN"] = os.environ[token_env_name]
        notes.append(f"TOWER_ACCESS_TOKEN=env:{token_env_name}")
    elif os.environ.get("TOWER_ACCESS_TOKEN"):
        env["TOWER_ACCESS_TOKEN"] = os.environ["TOWER_ACCESS_TOKEN"]
        notes.append("TOWER_ACCESS_TOKEN=env:TOWER_ACCESS_TOKEN")
    else:
        notes.append("TOWER_ACCESS_TOKEN=not_set")

    tower_workspace_id = getattr(args, "tower_workspace_id", None)
    if tower_workspace_id:
        env["TOWER_WORKSPACE_ID"] = str(tower_workspace_id)
        notes.append("TOWER_WORKSPACE_ID=set")

    tower_api_endpoint = getattr(args, "tower_api_endpoint", None)
    if tower_api_endpoint:
        env["TOWER_API_ENDPOINT"] = str(tower_api_endpoint)
        notes.append("TOWER_API_ENDPOINT=set")

    return notes


def state_path(args: argparse.Namespace) -> Path:
    return resolve_in_pipeline(args.state_file)


def load_state(args: argparse.Namespace) -> Dict[str, object]:
    path = state_path(args)
    if not path.exists():
        return {"created_at": now_iso(), "runs": {}}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Could not parse state file {path}: {exc}") from exc


def save_state(state: Dict[str, object], args: argparse.Namespace) -> None:
    path = state_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_manifest(state, args)


def history_info_from_command(command: str, args: argparse.Namespace) -> Optional[str]:
    del args
    # Match longer names first so barcode_filter_knee2 is not classified as
    # barcode_filter_knee. This intentionally follows the newest Nextflow
    # command for each run name, even if it was launched with a different local
    # state/config directory.
    for spec in sorted(RUN_SPECS, key=lambda item: len(item.info), reverse=True):
        if spec.info in command:
            return spec.info
    return None


def parse_nextflow_history_line(line: str, args: argparse.Namespace) -> Optional[Dict[str, str]]:
    fields = line.rstrip("\n").split("\t")
    command = fields[-1].strip() if fields else ""
    if "nextflow run main.nf" not in command or "-preview" in command:
        return None

    info = history_info_from_command(command, args)
    if info is None:
        return None

    session = ""
    if len(fields) >= 6:
        session = fields[-2].strip()
    if not NEXTFLOW_SESSION_RE.fullmatch(session):
        matches = NEXTFLOW_SESSION_RE.findall(line)
        session = matches[-1] if matches else ""
    if not session:
        return None

    status = fields[3].strip() if len(fields) >= 4 else ""
    return {
        "info": info,
        "started": fields[0].strip() if fields else "",
        "session": session,
        "history_status": status,
        "command": command,
    }


def read_nextflow_history(args: argparse.Namespace) -> List[Dict[str, str]]:
    env = os.environ.copy()
    env["NXF_VER"] = args.nxf_ver
    completed = subprocess.run(
        ["nextflow", "log", "-f", "session,status,command"],
        cwd=str(SCRIPT_DIR),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(
            "Could not read Nextflow history with `nextflow log`: "
            + completed.stderr.strip()
        )
    records = []
    for line in completed.stdout.splitlines():
        parsed = parse_nextflow_history_line(line, args)
        if parsed:
            records.append(parsed)
    records.sort(key=lambda record: record.get("started", ""))
    return records


def recover_history(state: Dict[str, object], args: argparse.Namespace) -> int:
    ensure_state_records(state, args)
    runs = state.get("runs", {})
    if not isinstance(runs, dict):
        return 0

    recovered = 0
    for history in read_nextflow_history(args):
        record = runs.get(history["info"])
        if not isinstance(record, dict):
            continue
        pid = record.get("pid")
        if isinstance(pid, int) and pid_is_alive(pid):
            continue
        history_status = history["history_status"]
        if history_status == "OK":
            record["last_ok_resume_id"] = history["session"]
        record["resume_id"] = history["session"]
        record["history_status"] = history_status
        record["history_command"] = history["command"]
        record["history_started_at"] = history.get("started", "")
        record["history_recovered_at"] = now_iso()
        mapped_status = NEXTFLOW_STATUS_TO_STATE.get(history_status)
        if mapped_status:
            record["status"] = mapped_status
            if mapped_status != RUNNING:
                record["pid"] = None
        recovered += 1

    save_state(state, args)
    if recovered:
        append_event(args, "history_recovered", extra={"recovered_records": recovered})
    return recovered


def ensure_state_records(state: Dict[str, object], args: argparse.Namespace) -> None:
    runs = state.setdefault("runs", {})
    assert isinstance(runs, dict)
    config_dir = resolve_in_pipeline(args.config_dir)
    log_dir = resolve_in_pipeline(args.log_dir)

    for spec in RUN_SPECS:
        record = runs.setdefault(spec.info, {})
        assert isinstance(record, dict)
        record.setdefault("status", PENDING)
        record.setdefault("attempts", 0)
        record["info"] = spec.info
        record["qc_barcode_filter"] = spec.qc_barcode_filter
        record["qc_min_genes_per_cell"] = spec.qc_min_genes_per_cell
        record["config_path"] = str(config_path_for(spec, config_dir))
        record["log_path"] = str(log_path_for(spec, log_dir))
        record["event_log_path"] = str(event_log_path(args))
        record["manifest_path"] = str(manifest_path(args))
        record["state_path"] = str(state_path(args))
        record["outdir"] = outdir_for(spec, args)
        if record.get("status", PENDING) == PENDING:
            record["command"] = display_command(build_command(spec, args), args)
            record["used_resume"] = False
            record["forced_relaunch"] = False
        else:
            record.setdefault("command", display_command(build_command(spec, args), args))
            record.setdefault("used_resume", False)
            record.setdefault("forced_relaunch", False)


def update_statuses(
    state: Dict[str, object],
    args: argparse.Namespace,
    processes: Optional[Dict[str, subprocess.Popen]] = None,
) -> bool:
    changed = False
    runs = state.get("runs", {})
    if not isinstance(runs, dict):
        return False

    local_pids = discover_local_nextflow_processes(args)
    for spec in RUN_SPECS:
        record = runs.get(spec.info)
        if not isinstance(record, dict):
            continue

        local_process = local_pids.get(spec.info)
        if local_process:
            local_pid = local_process.get("pid")
            local_command = local_process.get("command")
            if record.get("status") != RUNNING or record.get("pid") != local_pid:
                record["status"] = RUNNING
                record["pid"] = local_pid
                changed = True
            if isinstance(local_command, str):
                record["live_command"] = local_command
                if " -resume" in local_command:
                    record["used_resume"] = True
            continue

        if record.get("status") != RUNNING:
            continue

        process = processes.get(spec.info) if processes else None
        if process is not None:
            return_code = process.poll()
            if return_code is None:
                continue

            record["returncode"] = return_code
            record["ended_at"] = now_iso()
            record["status"] = SUCCEEDED if return_code == 0 else INTERRUPTED
            record["pid"] = None
            append_event(args, "run_finished", spec, record)
            changed = True
            continue

        pid = record.get("pid")
        pid_value = int(pid) if isinstance(pid, int) else None
        if pid_is_alive(pid_value):
            continue

        inferred = infer_status_from_log(Path(str(record.get("log_path", ""))))

        if record.get("history_status") == "-" and inferred is None:
            record["status"] = RUNNING
            record["pid"] = None
            changed = True
            continue

        record["status"] = inferred or INTERRUPTED
        record["ended_at"] = now_iso()
        record["pid"] = None
        append_event(args, "run_status_recovered", spec, record)
        changed = True

    if changed:
        save_state(state, args)
    return changed


def launch_run(
    spec: RunSpec,
    state: Dict[str, object],
    args: argparse.Namespace,
    processes: Optional[Dict[str, subprocess.Popen]] = None,
    *,
    resume: bool = False,
    force: bool = False,
) -> str:
    ensure_state_records(state, args)
    runs = state["runs"]
    assert isinstance(runs, dict)
    record = runs[spec.info]
    assert isinstance(record, dict)

    if record.get("status") == RUNNING:
        if pid_is_alive(record.get("pid")):
            return f"{spec.info} is already running with pid {record.get('pid')}"
        if record.get("history_status") == "-" and not force:
            return f"{spec.info} appears running in Nextflow history; not launching again"
    if record.get("status") == SUCCEEDED and not force:
        return f"{spec.info} already succeeded; not launching again"

    generate_configs(args)
    log_dir = resolve_in_pipeline(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_path_for(spec, log_dir)
    if force and resume and record.get("last_ok_resume_id"):
        previous_resume_id = record.get("last_ok_resume_id")
    else:
        previous_resume_id = record.get("resume_id")
    resume_id = args.resume_id or (str(previous_resume_id) if previous_resume_id else None)
    command = build_command(spec, args, resume=resume, resume_id=resume_id)
    command_text = display_command(command, args)

    attempts = int(record.get("attempts", 0)) + 1
    record.update(
        {
            "attempts": attempts,
            "started_at": now_iso(),
            "ended_at": None,
            "returncode": None,
            "status": DRY_RUN if args.dry_run else RUNNING,
            "pid": None,
            "command": command_text,
            "used_resume": resume,
            "forced_relaunch": force,
            "log_path": str(log_path),
            "event_log_path": str(event_log_path(args)),
            "manifest_path": str(manifest_path(args)),
            "state_path": str(state_path(args)),
            "outdir": outdir_for(spec, args),
            "resume_id": resume_id if resume else None,
        }
    )

    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write("\n" + "=" * 88 + "\n")
        log_handle.write(f"Attempt: {attempts}\n")
        log_handle.write(f"Started: {record['started_at']}\n")
        log_handle.write(f"info: {spec.info}\n")
        log_handle.write(f"QC_barcode_filter: {spec.qc_barcode_filter}\n")
        log_handle.write(f"QC_min_genes_per_cell: {spec.min_genes_display}\n")
        log_handle.write(f"Resume mode: {resume}\n")
        log_handle.write(f"Forced relaunch: {force}\n")
        if resume and resume_id:
            log_handle.write(f"Resume id: {resume_id}\n")
        log_handle.write(f"Command: {command_text}\n")
        log_handle.write(f"State file: {state_path(args)}\n")
        log_handle.write(f"Event log: {event_log_path(args)}\n")
        log_handle.write(f"Manifest: {manifest_path(args)}\n")
        log_handle.flush()

        if not args.dry_run:
            env = os.environ.copy()
            env["NXF_VER"] = args.nxf_ver
            try:
                tower_notes = configure_tower_env(env, args)
            except RuntimeError as exc:
                log_handle.write(f"Launch failed: Tower credentials: {exc}\n")
                record["status"] = INTERRUPTED
                record["ended_at"] = now_iso()
                append_event(args, "launch_failed", spec, record, {"error": str(exc)})
                save_state(state, args)
                return f"launch failed for {spec.info}: {exc}"
            log_handle.write(f"Tower credentials: {', '.join(tower_notes)}\n")
            log_handle.flush()
            try:
                process = subprocess.Popen(
                    command,
                    cwd=str(SCRIPT_DIR),
                    env=env,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except OSError as exc:
                log_handle.write(f"Launch failed: {exc}\n")
                record["status"] = INTERRUPTED
                record["ended_at"] = now_iso()
                append_event(args, "launch_failed", spec, record, {"error": str(exc)})
                save_state(state, args)
                return f"launch failed for {spec.info}: {exc}"
            record["pid"] = process.pid
            if processes is not None:
                processes[spec.info] = process

    if args.dry_run:
        event = "dry_run_recorded"
    elif resume:
        event = "run_resumed"
    elif force:
        event = "run_force_relaunched"
    else:
        event = "run_launched"
    append_event(args, event, spec, record)
    save_state(state, args)
    if args.dry_run:
        return f"dry-run recorded for {spec.info}"
    if resume:
        return f"resumed {spec.info} with pid {record['pid']}"
    if force:
        return f"force relaunched {spec.info} with pid {record['pid']}"
    return f"launched {spec.info} with pid {record['pid']}"


def launch_many(
    specs: Iterable[RunSpec],
    state: Dict[str, object],
    args: argparse.Namespace,
    processes: Optional[Dict[str, subprocess.Popen]] = None,
    *,
    resume_interrupted: bool = False,
    resume_all: bool = False,
    force: bool = False,
) -> List[str]:
    messages = []
    for spec in specs:
        runs = state.get("runs", {})
        record = runs.get(spec.info, {}) if isinstance(runs, dict) else {}
        status = record.get("status", PENDING) if isinstance(record, dict) else PENDING
        use_resume = resume_all or (resume_interrupted and status in RESUMABLE_STATUSES)
        messages.append(
            launch_run(spec, state, args, processes, resume=use_resume, force=force)
        )
    return messages


def record_used_resume(record: Dict[str, object]) -> bool:
    if record.get("used_resume", False):
        return True
    for key in ("command", "history_command", "live_command"):
        value = record.get(key)
        if isinstance(value, str) and " -resume" in value:
            return True
    return False


def status_rows(state: Dict[str, object]) -> List[Dict[str, object]]:
    runs = state.get("runs", {})
    if not isinstance(runs, dict):
        return []
    rows = []
    for spec in RUN_SPECS:
        record = runs.get(spec.info, {})
        if not isinstance(record, dict):
            record = {}
        rows.append(
            {
                "info": spec.info,
                "filter": spec.qc_barcode_filter,
                "min_genes": spec.min_genes_display,
                "status": record.get("status", PENDING),
                "pid": record.get("pid") or "",
                "attempts": record.get("attempts", 0),
                "used_resume": record_used_resume(record),
                "forced_relaunch": record.get("forced_relaunch", False),
                "outdir": record.get("outdir", ""),
                "log_path": record.get("log_path", ""),
            }
        )
    return rows


def print_status(state: Dict[str, object]) -> None:
    rows = status_rows(state)
    print(f"{'status':<12} {'info':<32} {'filter':<8} {'min_genes':<9} {'pid':<8} {'resume':<6} attempts")
    print("-" * 90)
    for row in rows:
        print(
            f"{str(row['status']):<12} "
            f"{str(row['info']):<32} "
            f"{str(row['filter']):<8} "
            f"{str(row['min_genes']):<9} "
            f"{str(row['pid']):<8} "
            f"{str(row['used_resume']):<6} "
            f"{row['attempts']}"
        )


def print_commands(args: argparse.Namespace) -> None:
    generate_configs(args)
    for spec in RUN_SPECS:
        print(f"[{spec.info}]")
        print(display_command(build_command(spec, args), args))
        print(f"config: {config_path_for(spec, resolve_in_pipeline(args.config_dir))}")
        print(f"outdir: {outdir_for(spec, args)}")
        print()


def addnstr_safe(stdscr: curses.window, y: int, x: int, text: str, width: int, attr: int = 0) -> None:
    height, screen_width = stdscr.getmaxyx()
    if y < 0 or y >= height or x >= screen_width:
        return
    clipped_width = max(0, min(width, screen_width - x))
    if clipped_width == 0:
        return
    stdscr.addnstr(y, x, text.ljust(clipped_width), clipped_width, attr)


class SweepTui:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.state = load_state(args)
        ensure_state_records(self.state, args)
        save_state(self.state, args)
        self.processes: Dict[str, subprocess.Popen] = {}
        self.selected = 0
        self.interrupted_only = False
        self.message = "s fresh-starts selected; f force-relaunches completed; r resumes interrupted"
        self.pending_start_all = False
        self.last_status_save = 0.0

    def visible_specs(self) -> List[RunSpec]:
        if not self.interrupted_only:
            return list(RUN_SPECS)
        runs = self.state.get("runs", {})
        if not isinstance(runs, dict):
            return []
        visible = []
        for spec in RUN_SPECS:
            record = runs.get(spec.info, {})
            if isinstance(record, dict) and record.get("status") in RESUMABLE_STATUSES:
                visible.append(spec)
        return visible

    def selected_spec(self) -> Optional[RunSpec]:
        visible = self.visible_specs()
        if not visible:
            return None
        self.selected = max(0, min(self.selected, len(visible) - 1))
        return visible[self.selected]

    def refresh_statuses(self) -> None:
        if time.time() - self.last_status_save < 1.0:
            return
        update_statuses(self.state, self.args, self.processes)
        self.last_status_save = time.time()

    def start_selected(self) -> None:
        spec = self.selected_spec()
        if spec is None:
            self.message = "no run selected"
            return
        self.message = launch_run(spec, self.state, self.args, self.processes)

    def force_relaunch_selected(self) -> None:
        spec = self.selected_spec()
        if spec is None:
            self.message = "no run selected"
            return
        self.message = launch_run(spec, self.state, self.args, self.processes, force=True)

    def resume_selected_or_show_interrupted(self) -> None:
        spec = self.selected_spec()
        runs = self.state.get("runs", {})
        record = runs.get(spec.info, {}) if isinstance(runs, dict) and spec else {}
        if spec is not None and isinstance(record, dict) and record.get("status") in RESUMABLE_STATUSES:
            self.message = launch_run(spec, self.state, self.args, self.processes, resume=True)
            return

        interrupted = [
            run_spec
            for run_spec in RUN_SPECS
            if isinstance(runs, dict)
            and isinstance(runs.get(run_spec.info, {}), dict)
            and runs[run_spec.info].get("status") in RESUMABLE_STATUSES
        ]
        if not interrupted:
            self.message = "no interrupted runs found"
            self.interrupted_only = False
            return
        self.interrupted_only = True
        self.selected = 0
        self.message = "showing interrupted runs; press r again to resume the selected one"

    def start_all_confirmed(self) -> None:
        if not self.pending_start_all:
            self.pending_start_all = True
            self.message = "press a again within 5 seconds to start pending and resume interrupted runs"
            self.start_all_requested_at = time.time()
            return

        if time.time() - getattr(self, "start_all_requested_at", 0.0) > 5.0:
            self.pending_start_all = False
            self.message = "start-all confirmation expired; press a again"
            return

        runs = self.state.get("runs", {})
        specs = []
        for spec in RUN_SPECS:
            if spec.info in set(self.args.exclude or []):
                continue
            record = runs.get(spec.info, {}) if isinstance(runs, dict) else {}
            status = record.get("status", PENDING) if isinstance(record, dict) else PENDING
            if status not in {RUNNING, SUCCEEDED}:
                specs.append(spec)

        messages = launch_many(specs, self.state, self.args, self.processes, resume_interrupted=True)
        self.pending_start_all = False
        self.message = "; ".join(messages[-2:]) if messages else "nothing to launch"

    def draw(self, stdscr: curses.window) -> None:
        stdscr.erase()
        height, width = stdscr.getmaxyx()
        title = "QC sweep Nextflow TUI"
        if self.args.dry_run:
            title += " [dry-run]"
        if self.interrupted_only:
            title += " [interrupted only]"
        addnstr_safe(stdscr, 0, 0, title, width, curses.A_BOLD)
        addnstr_safe(
            stdscr,
            1,
            0,
            "q quit | up/down select | g regenerate configs | s fresh start | f force relaunch | a start/resume all | r show/resume interrupted",
            width,
        )
        addnstr_safe(stdscr, 2, 0, self.message, width, curses.A_DIM)

        header = f"{'status':<12} {'info':<32} {'filter':<8} {'min_genes':<9} {'pid':<8} {'resume':<6} tries"
        addnstr_safe(stdscr, 4, 0, header, width, curses.A_UNDERLINE)

        visible = self.visible_specs()
        if not visible:
            addnstr_safe(stdscr, 6, 0, "No runs match this view.", width)
        rows_by_info = {row["info"]: row for row in status_rows(self.state)}
        max_rows = max(0, height - 13)
        for idx, spec in enumerate(visible[:max_rows]):
            row = rows_by_info[spec.info]
            attr = curses.A_REVERSE if idx == self.selected else 0
            if row["status"] == SUCCEEDED:
                attr |= curses.A_BOLD
            text = (
                f"{str(row['status']):<12} "
                f"{spec.info:<32} "
                f"{spec.qc_barcode_filter:<8} "
                f"{spec.min_genes_display:<9} "
                f"{str(row['pid']):<8} "
                f"{str(row['used_resume']):<6} "
                f"{row['attempts']}"
            )
            addnstr_safe(stdscr, 5 + idx, 0, text, width, attr)

        selected = self.selected_spec()
        panel_y = max(7, min(height - 7, 6 + len(visible[:max_rows]) + 1))
        addnstr_safe(stdscr, panel_y, 0, "-" * width, width)

        if selected is None:
            stdscr.refresh()
            return

        runs = self.state.get("runs", {})
        record = runs.get(selected.info, {}) if isinstance(runs, dict) else {}
        if not isinstance(record, dict):
            record = {}
        addnstr_safe(stdscr, panel_y + 1, 0, f"outdir: {record.get('outdir', '')}", width)
        addnstr_safe(stdscr, panel_y + 2, 0, f"config: {record.get('config_path', '')}", width)
        addnstr_safe(stdscr, panel_y + 3, 0, f"log:    {record.get('log_path', '')}", width)

        log_path = Path(str(record.get("log_path", "")))
        available_tail_lines = max(0, height - panel_y - 5)
        for offset, line in enumerate(tail_lines(log_path, available_tail_lines)):
            addnstr_safe(stdscr, panel_y + 4 + offset, 0, line, width)
        stdscr.refresh()

    def handle_key(self, key: int) -> bool:
        visible = self.visible_specs()
        if key in (ord("q"), ord("Q")):
            return False
        if key in (curses.KEY_UP, ord("k")):
            self.selected = max(0, self.selected - 1)
        elif key in (curses.KEY_DOWN, ord("j")):
            self.selected = min(max(0, len(visible) - 1), self.selected + 1)
        elif key == ord("g"):
            written = generate_configs(self.args)
            ensure_state_records(self.state, self.args)
            save_state(self.state, self.args)
            self.message = f"generated {len(written)} configs in {resolve_in_pipeline(self.args.config_dir)}"
        elif key == ord("s"):
            self.start_selected()
        elif key == ord("f"):
            self.force_relaunch_selected()
        elif key == ord("a"):
            self.start_all_confirmed()
        elif key == ord("r"):
            self.resume_selected_or_show_interrupted()
        elif key == ord("v"):
            self.interrupted_only = not self.interrupted_only
            self.selected = 0
            self.message = "toggled interrupted-only view"
        else:
            self.pending_start_all = False
        return True

    def run(self, stdscr: curses.window) -> None:
        curses.curs_set(0)
        stdscr.nodelay(True)
        stdscr.timeout(300)

        while True:
            self.refresh_statuses()
            self.draw(stdscr)
            key = stdscr.getch()
            if key == -1:
                continue
            if not self.handle_key(key):
                break


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate, submit, and monitor six CRISPR pipeline QC sweep runs."
    )
    parser.add_argument("--base-config", default=DEFAULT_BASE_CONFIG)
    parser.add_argument("--config-dir", default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--log-dir", default=DEFAULT_LOG_DIR)
    parser.add_argument("--state-file", default=DEFAULT_STATE_FILE)
    parser.add_argument("--event-log", default=DEFAULT_EVENT_LOG)
    parser.add_argument("--manifest-file", default=DEFAULT_MANIFEST_FILE)
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--outdir-template", default=DEFAULT_OUTDIR_TEMPLATE)
    parser.add_argument(
        "--resume-id",
        default=DEFAULT_RESUME_ID,
        help="Optional Nextflow resume id used only for resume attempts.",
    )
    parser.add_argument("--nxf-ver", default=DEFAULT_NXF_VER)
    parser.add_argument("--profile", default="google")
    tower_group = parser.add_mutually_exclusive_group()
    tower_group.add_argument(
        "--with-tower",
        dest="with_tower",
        action="store_true",
        default=DEFAULT_WITH_TOWER,
        help="include Nextflow --with_tower in launched commands (default)",
    )
    tower_group.add_argument(
        "--without-tower",
        dest="with_tower",
        action="store_false",
        help="omit Nextflow --with_tower from launched commands",
    )
    parser.add_argument(
        "--tower-access-token-env",
        default=DEFAULT_TOWER_ACCESS_TOKEN_ENV,
        help=(
            "environment variable to copy into TOWER_ACCESS_TOKEN for launched "
            "Nextflow processes"
        ),
    )
    parser.add_argument(
        "--tower-access-token-file",
        help=(
            "file containing the Tower/Seqera access token; value is passed as "
            "TOWER_ACCESS_TOKEN and is not written to logs"
        ),
    )
    parser.add_argument(
        "--tower-workspace-id",
        help="optional Tower/Seqera workspace id passed as TOWER_WORKSPACE_ID",
    )
    parser.add_argument(
        "--tower-api-endpoint",
        help="optional Tower/Seqera API endpoint passed as TOWER_API_ENDPOINT",
    )
    parser.add_argument("--dry-run", action="store_true", help="print/record commands without starting Nextflow")
    parser.add_argument("--generate-only", action="store_true", help="write the six config files and exit")
    parser.add_argument("--status", action="store_true", help="print saved status and exit")
    parser.add_argument("--commands", action="store_true", help="print all generated commands and exit")
    parser.add_argument(
        "--launch",
        choices=[spec.info for spec in RUN_SPECS],
        help="launch one run without opening the TUI",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="with --launch, add Nextflow -resume; TUI fresh starts are unaffected",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="with --launch, fresh relaunch even if saved status is succeeded; still refuses active running PIDs",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        choices=[spec.info for spec in RUN_SPECS],
        help="exclude a run from --launch-all; can be used more than once",
    )
    parser.add_argument(
        "--launch-all",
        action="store_true",
        help="launch pending runs and resume interrupted runs without the TUI",
    )
    parser.add_argument(
        "--recover-history-only",
        action="store_true",
        help="read Nextflow history into the TUI state and exit",
    )
    parser.add_argument("--no-tui", action="store_true", help="do not open curses TUI")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    generated = generate_configs(args)
    append_event(
        args,
        "configs_generated",
        extra={
            "config_count": len(generated),
            "config_dir": str(resolve_in_pipeline(args.config_dir)),
        },
    )
    state = load_state(args)
    ensure_state_records(state, args)
    recovered = recover_history(state, args)
    update_statuses(state, args)
    save_state(state, args)

    if args.recover_history_only:
        print(f"Recovered {recovered} Nextflow history records into {state_path(args)}")
        return 0

    if args.commands or args.dry_run:
        print_commands(args)
        return 0

    if args.generate_only:
        print(f"Generated {len(RUN_SPECS)} configs in {resolve_in_pipeline(args.config_dir)}")
        return 0

    if args.status:
        print_status(state)
        return 0

    if args.launch:
        message = launch_run(spec_by_info(args.launch), state, args, resume=args.resume, force=args.force)
        print(message)
        return 0

    if args.launch_all:
        runs = state.get("runs", {})
        excluded = set(args.exclude or [])
        specs = []
        for spec in RUN_SPECS:
            if spec.info in excluded:
                continue
            record = runs.get(spec.info, {}) if isinstance(runs, dict) else {}
            status = record.get("status", PENDING) if isinstance(record, dict) else PENDING
            if args.force:
                if status != RUNNING:
                    specs.append(spec)
            elif status not in {RUNNING, SUCCEEDED}:
                specs.append(spec)
        for message in launch_many(
            specs,
            state,
            args,
            resume_interrupted=True,
            resume_all=args.resume,
            force=args.force,
        ):
            print(message)
        return 0

    if args.no_tui:
        print_status(state)
        return 0

    if not sys.stdin.isatty() or not sys.stdout.isatty():
        print_status(state)
        print("\nNo TTY detected; rerun without redirection to open the TUI.")
        return 0

    curses.wrapper(lambda stdscr: SweepTui(args).run(stdscr))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
