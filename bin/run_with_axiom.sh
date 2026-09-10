#!/usr/bin/env bash
# Run any Nextflow command with fail-open Axiom process telemetry.
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_NAME="${AXIOM_RUN_NAME:-crispr_$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_ID="${AXIOM_RUN_ID:-$(python -c 'import uuid; print(uuid.uuid4())')}"
STATE_DIR="${AXIOM_STATE_DIR:-$PWD/.axiom_telemetry/$RUN_ID}"
OUTDIR="${AXIOM_OUTDIR:?Set AXIOM_OUTDIR to the same output directory passed to Nextflow}"
TRACE_FILE="$STATE_DIR/nextflow_trace.tsv"
HOOK_FILE="$STATE_DIR/nextflow_hooks.jsonl"
DAG_FILE="$STATE_DIR/nextflow_dag.dot"

mkdir -p "$STATE_DIR"
[[ -n "${AXIOM_IGVF:-}" ]] || { echo "AXIOM_IGVF is not set" >&2; exit 2; }
[[ $# -ge 2 && "$(basename "$1")" == "nextflow" ]] || {
  echo "usage: AXIOM_OUTDIR=... run_with_axiom.sh nextflow run <pipeline> [options]" >&2
  exit 2
}
HAS_RUN=false
for argument in "${@:2}"; do
  [[ "$argument" == "run" ]] && HAS_RUN=true
done
[[ "$HAS_RUN" == "true" ]] || {
  echo "run_with_axiom.sh requires a Nextflow run command" >&2
  exit 2
}

NEXTFLOW_BIN="$1"
shift
DASHBOARD_FLAG="--create-dashboard"
[[ "${AXIOM_CREATE_DASHBOARD:-true}" == "true" ]] || DASHBOARD_FLAG="--no-create-dashboard"
INGEST_URL="${AXIOM_INGEST_URL:-}"
[[ -n "$INGEST_URL" ]] || INGEST_URL='https://us-east-1.aws.edge.axiom.co/v1/ingest/{dataset}'

exec python "$PIPELINE_DIR/bin/axiom_telemetry.py" run \
  --dataset "${AXIOM_DATASET:-crispr-pipeline}" \
  --token-env "${AXIOM_TOKEN_ENV:-AXIOM_IGVF}" \
  --ingest-url "$INGEST_URL" \
  --api-url "${AXIOM_API_URL:-https://api.axiom.co}" \
  --max-bytes "${AXIOM_MAX_BYTES:-20000000}" \
  --run-id "$RUN_ID" --run-name "$RUN_NAME" \
  --pipeline-revision "$(git -C "$PIPELINE_DIR" rev-parse --short=12 HEAD 2>/dev/null || echo unknown)" \
  --state-dir "$STATE_DIR" --trace "$TRACE_FILE" --dag "$DAG_FILE" --hook-events "$HOOK_FILE" --outdir "$OUTDIR" \
  "$DASHBOARD_FLAG" \
  -- "$NEXTFLOW_BIN" -c "$PIPELINE_DIR/conf/axiom.config" "$@" \
  -with-trace "$TRACE_FILE" -with-dag "$DAG_FILE" \
  --AXIOM_run_id "$RUN_ID" --AXIOM_run_name "$RUN_NAME" --AXIOM_event_file "$HOOK_FILE"
