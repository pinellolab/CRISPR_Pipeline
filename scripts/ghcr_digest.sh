#!/usr/bin/env bash
# Print the manifest digest GHCR serves for a PerTurbo tag.
#
# ghcr.io requires a bearer token even for a public repository, but the token
# endpoint hands one out anonymously for pull scope.  The digest is whatever the
# registry returns in docker-content-digest for a HEAD on the manifest -- which
# is exactly the string nextflow.config should pin.
#
# Usage:  scripts/ghcr_digest.sh v2.0.0rc11
#         scripts/ghcr_digest.sh v2.0.0rc11 --json   # tag, digest, media type
set -Eeuo pipefail

REPOSITORY="${GHCR_REPOSITORY:-pinellolab/perturbo}"
REGISTRY="https://ghcr.io"

usage() {
    cat >&2 <<USAGE
usage: ${0##*/} TAG [--json]

Prints the docker-content-digest ghcr.io/${REPOSITORY} serves for TAG.
Override the repository with GHCR_REPOSITORY=owner/name.
USAGE
}

TAG=""
AS_JSON=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --json) AS_JSON=true ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "unknown option: $1" >&2; usage; exit 2 ;;
        *)
            if [[ -n "$TAG" ]]; then echo "unexpected argument: $1" >&2; usage; exit 2; fi
            TAG="$1"
            ;;
    esac
    shift
done
[[ -n "$TAG" ]] || { usage; exit 2; }

command -v curl >/dev/null || { echo "curl is required" >&2; exit 1; }

# Anonymous pull token.  No jq dependency: the payload is a flat JSON object.
token_json="$(curl -fsSL "${REGISTRY}/token?scope=repository:${REPOSITORY}:pull")" || {
    echo "failed to get an anonymous pull token for ${REPOSITORY}" >&2
    exit 1
}
TOKEN="$(printf '%s' "$token_json" | sed -n 's/.*"token"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p')"
[[ -n "$TOKEN" ]] || { echo "could not parse a token out of: $token_json" >&2; exit 1; }

# Accept every manifest form the tag might be: an OCI image index or a Docker
# manifest list for a multi-arch tag, a plain manifest for a single-arch one.
# Without these the registry answers with a schema-1 manifest whose digest is
# not the digest anyone else sees.
headers="$(curl -fsSI \
    -H "Authorization: Bearer ${TOKEN}" \
    -H "Accept: application/vnd.oci.image.index.v1+json" \
    -H "Accept: application/vnd.oci.image.manifest.v1+json" \
    -H "Accept: application/vnd.docker.distribution.manifest.v2+json" \
    -H "Accept: application/vnd.docker.distribution.manifest.list.v2+json" \
    "${REGISTRY}/v2/${REPOSITORY}/manifests/${TAG}")" || {
    echo "HEAD on ${REGISTRY}/v2/${REPOSITORY}/manifests/${TAG} failed (tag missing?)" >&2
    exit 1
}

# Header names are case-insensitive; HTTP/2 lowercases them, HTTP/1.1 may not.
DIGEST="$(printf '%s' "$headers" \
    | tr -d '\r' \
    | awk 'BEGIN{IGNORECASE=1} /^docker-content-digest:/ {print $2}' \
    | tail -1)"
MEDIA="$(printf '%s' "$headers" \
    | tr -d '\r' \
    | awk 'BEGIN{IGNORECASE=1} /^content-type:/ {print $2}' \
    | tail -1)"

if [[ ! "$DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]]; then
    echo "no usable docker-content-digest in the response for ${TAG}" >&2
    printf '%s\n' "$headers" >&2
    exit 1
fi

if [[ "$AS_JSON" == true ]]; then
    printf '{"repository":"%s","tag":"%s","digest":"%s","mediaType":"%s"}\n' \
        "$REPOSITORY" "$TAG" "$DIGEST" "$MEDIA"
else
    printf '%s\n' "$DIGEST"
fi
