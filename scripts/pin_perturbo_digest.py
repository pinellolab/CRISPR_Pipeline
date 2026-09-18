#!/usr/bin/env python3
"""Repin the PerTurbo container image across the pipeline's Nextflow configs.

The pipeline pins PerTurbo by immutable digest in ``nextflow.config``
(``params.containers.perturbo``).  Site configs may carry their own override --
another digest, a tag, or a local ``.sif`` built from that tag.  This script
rewrites all of them at once so a release bump cannot leave one behind, and
prints a unified diff of every line it touched.

It also appends a bullet to the ``### Inference`` section of ``CHANGELOG.md``
recording the new pin, with a placeholder sentence for a human to complete.

Usage:
    scripts/pin_perturbo_digest.py --tag v2.0.0rc11 --digest sha256:<64 hex>
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from pathlib import Path

IMAGE = "ghcr.io/pinellolab/perturbo"

# A digest must be the full, unambiguous thing: sha256: plus exactly 64 hex
# characters.  A short digest is not a pin -- the registry will not resolve it
# and Nextflow will not error until the task runs.
DIGEST_RE = re.compile(r"\Asha256:[0-9a-f]{64}\Z")

# Tag looks like v2.0.0rc11 / v2.0.0 / 2.0.0rc11.
TAG_RE = re.compile(r"\Av?\d+\.\d+\.\d+(?:rc\d+)?\Z")

# Config files that may carry a PerTurbo image reference.  Order is only for a
# stable report; every match is rewritten.
CONFIG_GLOBS = (
    "nextflow.config",
    "nextflow_*.config",
    "conf/*.config",
    "conf/*.conf",
)

IMAGE_RE = re.escape(IMAGE)

# 1. digest pin: ghcr.io/pinellolab/perturbo@sha256:<64 hex>
SUB_DIGEST = re.compile(rf"{IMAGE_RE}@sha256:[0-9a-fA-F]{{64}}")
# 2. tag pin: ghcr.io/pinellolab/perturbo:<tag>
SUB_TAG = re.compile(rf"{IMAGE_RE}:[^\s'\"]+")
# 3. a versioned local image file: perturbo_v2.0.0rc10.sif, perturbo-2.0.0rc9.sif
#    A bare perturbo.sif is deliberately NOT touched: reproduction bundles use
#    that fixed name and renaming it would break their launchers.
SUB_SIF = re.compile(r"perturbo[_-]v?\d+\.\d+\.\d+(?:rc\d+)?\.(sif|img)")
# 4. the human-readable comment that names the release above the pin,
#    e.g. "// perturbo v2.0 rc9"
SUB_COMMENT = re.compile(r"(//\s*perturbo\s+)v?\d[\w.]*\s*rc\d+", re.IGNORECASE)

CHANGELOG_HEADING = "### Inference"


def rewrite_text(text: str, tag: str, digest: str) -> str:
    """Apply every PerTurbo pin substitution to a config file's text."""
    out = SUB_DIGEST.sub(f"{IMAGE}@{digest}", text)
    # Only after the digest form is gone, so ':' in "@sha256:..." is not a tag.
    out = SUB_TAG.sub(f"{IMAGE}:{tag}", out)
    out = SUB_SIF.sub(lambda m: f"perturbo_{tag}.{m.group(1)}", out)
    out = SUB_COMMENT.sub(rf"\g<1>{tag}", out)
    return out


def iter_config_paths(repo: Path):
    seen = set()
    for pattern in CONFIG_GLOBS:
        for path in sorted(repo.glob(pattern)):
            if path.is_file() and path not in seen:
                seen.add(path)
                yield path


def changelog_bullet(tag: str, digest: str) -> str:
    return (
        f"- Pin PerTurbo to `{tag}` (`{digest}`). "
        "TODO: what changed in this PerTurbo release, and what was validated "
        "against it before the pin moved."
    )


def update_changelog(path: Path, tag: str, digest: str) -> tuple[str, str] | None:
    """Append the pin bullet to the ### Inference section.

    Returns (before, after) when the file changes, else None.
    """
    before = path.read_text()
    bullet = changelog_bullet(tag, digest)
    if f"Pin PerTurbo to `{tag}` (`{digest}`)" in before:
        return None  # already recorded; do not duplicate on a rerun

    lines = before.splitlines(keepends=True)
    try:
        start = next(
            i for i, line in enumerate(lines) if line.strip() == CHANGELOG_HEADING
        )
    except StopIteration:
        raise SystemExit(
            f"error: no '{CHANGELOG_HEADING}' heading in {path}; refusing to guess "
            "where the bullet belongs"
        )

    # The section runs to the next heading of the same or higher level.
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if re.match(r"^#{1,3} ", lines[i]):
            end = i
            break

    # Insert after the last top-level bullet in the section, so the new pin sits
    # next to the pin it replaces rather than after the prose that follows.
    last_bullet = None
    for i in range(start + 1, end):
        if lines[i].startswith("- "):
            last_bullet = i
    insert_at = (last_bullet + 1) if last_bullet is not None else (start + 1)

    lines.insert(insert_at, "\n" + bullet + "\n")
    after = "".join(lines)
    return (before, after)


def print_diff(path: Path, before: str, after: str, repo: Path) -> None:
    rel = path.relative_to(repo)
    diff = difflib.unified_diff(
        before.splitlines(keepends=True),
        after.splitlines(keepends=True),
        fromfile=f"a/{rel}",
        tofile=f"b/{rel}",
    )
    sys.stdout.writelines(diff)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Repin the PerTurbo container across the pipeline's configs."
    )
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v2.0.0rc11")
    parser.add_argument(
        "--digest",
        required=True,
        help="Image digest: sha256: followed by exactly 64 hex characters",
    )
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Repository root (default: the checkout this script lives in)",
    )
    parser.add_argument(
        "--no-changelog",
        action="store_true",
        help="Rewrite the configs only; leave CHANGELOG.md alone",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the diff without writing anything",
    )
    args = parser.parse_args(argv)

    digest = args.digest.strip()
    if not DIGEST_RE.match(digest):
        parser.error(
            f"--digest must be 'sha256:' followed by exactly 64 lowercase hex "
            f"characters; got {args.digest!r}"
        )
    tag = args.tag.strip()
    if not TAG_RE.match(tag):
        parser.error(f"--tag does not look like a release tag: {args.tag!r}")

    repo = args.repo.resolve()
    if not (repo / "nextflow.config").is_file():
        parser.error(f"{repo} does not look like the pipeline repo (no nextflow.config)")

    changed = []
    seen_reference = False
    for path in iter_config_paths(repo):
        before = path.read_text()
        if SUB_DIGEST.search(before) or SUB_TAG.search(before) or SUB_SIF.search(before):
            seen_reference = True
        after = rewrite_text(before, tag, digest)
        if after != before:
            changed.append(path)
            print_diff(path, before, after, repo)
            if not args.dry_run:
                path.write_text(after)

    if not seen_reference:
        print(
            "warning: no PerTurbo image reference found in any config; nothing repinned",
            file=sys.stderr,
        )
    elif not changed:
        print(
            f"note: every PerTurbo config reference already names {tag} ({digest})",
            file=sys.stderr,
        )

    if not args.no_changelog:
        changelog = repo / "CHANGELOG.md"
        if not changelog.is_file():
            print(f"warning: {changelog} not found; skipping", file=sys.stderr)
        else:
            result = update_changelog(changelog, tag, digest)
            if result is None:
                print(
                    f"note: CHANGELOG.md already records {tag} ({digest}); not duplicating",
                    file=sys.stderr,
                )
            else:
                before, after = result
                changed.append(changelog)
                print_diff(changelog, before, after, repo)
                if not args.dry_run:
                    changelog.write_text(after)

    verb = "would change" if args.dry_run else "changed"
    print(
        f"\n{verb} {len(changed)} file(s): "
        + ", ".join(str(p.relative_to(repo)) for p in changed),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
