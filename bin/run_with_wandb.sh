#!/usr/bin/env bash
# Run Nextflow with a fail-open, single-HTML W&B telemetry sidecar.
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_NAME="${WANDB_RUN_NAME:-crispr_$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_ID="${WANDB_SOURCE_RUN_ID:-$(python -c 'import uuid; print(uuid.uuid4())')}"
STATE_DIR="${WANDB_STATE_DIR:-$PWD/.wandb_telemetry/$RUN_ID}"
OUTDIR="${WANDB_OUTDIR:?Set WANDB_OUTDIR to the same output directory passed to Nextflow}"
TRACE_FILE="$STATE_DIR/nextflow_trace.tsv"
NEXTFLOW_LOG="$STATE_DIR/nextflow.log"
STATUS_FILE="$STATE_DIR/status.json"
DASHBOARD_HTML="$STATE_DIR/pipeline_execution.html"
PYTHON_BIN="${WANDB_PYTHON:-python}"

mkdir -p "$STATE_DIR"
[[ $# -ge 2 && "$(basename "$1")" == "nextflow" ]] || {
  echo "usage: WANDB_OUTDIR=... run_with_wandb.sh nextflow run <pipeline> [options]" >&2
  exit 2
}

NEXTFLOW_BIN="$1"
shift
HAS_RUN=false
for argument in "$@"; do
  [[ "$argument" == "run" ]] && HAS_RUN=true
done
[[ "$HAS_RUN" == "true" ]] || {
  echo "run_with_wandb.sh requires a Nextflow run command" >&2
  exit 2
}

"$PYTHON_BIN" "$PIPELINE_DIR/bin/wandb_html_monitor.py" \
  --project "${WANDB_PROJECT:-crispr-pipeline}" \
  --entity "${WANDB_ENTITY:-}" \
  --token-env "${WANDB_TOKEN_ENV:-WB_IGVF}" \
  --run-name "$RUN_NAME" --source-run-id "$RUN_ID" \
  --outdir "$OUTDIR" --trace "$TRACE_FILE" --nextflow-log "$NEXTFLOW_LOG" \
  --status-file "$STATUS_FILE" --dashboard-html "$DASHBOARD_HTML" \
  --poll-seconds "${WANDB_POLL_SECONDS:-30}" \
  --max-total-bytes "${WANDB_MAX_BYTES:-20000000}" &
MONITOR_PID=$!

write_status() {
  local status="$1"
  local exit_code="$2"
  printf '{"status":"%s","exit_code":%s}\n' "$status" "$exit_code" > "$STATUS_FILE"
}

finish_monitor() {
  local wrapper_exit=$?
  if [[ ! -f "$STATUS_FILE" ]]; then
    write_status interrupted "$wrapper_exit"
  fi
  wait "$MONITOR_PID" || echo "WARN: W&B HTML telemetry stopped; Nextflow exit code is unchanged" >&2
}
trap finish_monitor EXIT

set +e
"$NEXTFLOW_BIN" -log "$NEXTFLOW_LOG" "$@" -with-trace "$TRACE_FILE"
NEXTFLOW_EXIT=$?
set -e

RUN_STATUS=failed
[[ "$NEXTFLOW_EXIT" -eq 0 ]] && RUN_STATUS=completed
write_status "$RUN_STATUS" "$NEXTFLOW_EXIT"

wait "$MONITOR_PID" || true
trap - EXIT
exit "$NEXTFLOW_EXIT"
