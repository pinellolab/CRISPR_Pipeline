#!/usr/bin/env bash
# One-time offline installation of the versioned dashboard plugin archive.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION=0.1.2
ARCHIVE="$ROOT/plugins/nf-crispr-dashboard/dist/nf-crispr-dashboard-$VERSION.zip"
EXPECTED_SHA256=02a1eba42b8f30dc29952a603a36c485a089dfa5ea18d7104e415040ec902bbf

[[ -n "${NXF_PLUGINS_DIR:-}" && "$NXF_PLUGINS_DIR" = /* ]] || {
    echo 'Set NXF_PLUGINS_DIR to an absolute, writable shared plugin-cache path.' >&2
    exit 2
}
if command -v sha256sum >/dev/null 2>&1; then
    ACTUAL_SHA256="$(sha256sum "$ARCHIVE" | cut -d ' ' -f 1)"
else
    ACTUAL_SHA256="$(shasum -a 256 "$ARCHIVE" | cut -d ' ' -f 1)"
fi
[[ "$ACTUAL_SHA256" == "$EXPECTED_SHA256" ]] || {
    echo 'Dashboard plugin archive checksum mismatch; refusing installation.' >&2
    exit 1
}
TARGET="$NXF_PLUGINS_DIR/nf-crispr-dashboard-$VERSION"
if [[ -e "$TARGET" ]]; then
    [[ -f "$TARGET/.archive.sha256" && "$(<"$TARGET/.archive.sha256")" == "$EXPECTED_SHA256" ]] || {
        echo "Existing plugin at $TARGET differs; refusing to overwrite it." >&2
        exit 1
    }
    echo "Dashboard plugin already installed: $TARGET"
    exit 0
fi

mkdir -p "$NXF_PLUGINS_DIR"
STAGING="$(mktemp -d "$NXF_PLUGINS_DIR/.nf-crispr-dashboard.XXXXXX")"
trap 'rmdir "$STAGING" 2>/dev/null || true' EXIT
unzip -q "$ARCHIVE" -d "$STAGING"
printf '%s\n' "$EXPECTED_SHA256" > "$STAGING/.archive.sha256"
mv "$STAGING" "$TARGET"
trap - EXIT
echo "Installed dashboard plugin: $TARGET"
