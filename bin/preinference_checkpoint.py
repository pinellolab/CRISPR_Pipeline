#!/usr/bin/env python3
"""Publish an immutable MuData handoff immediately before inference."""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import mudata as mu


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", default="preinference_mudata.h5mu", type=Path)
    parser.add_argument("--manifest", default="preinference_manifest.json", type=Path)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()

    if args.output.resolve() == args.input.resolve():
        parser.error("checkpoint output must not overwrite the input")
    shutil.copyfile(args.input, args.output)
    mdata = mu.read_h5mu(args.output, backed="r")
    try:
        if not {"gene", "guide"}.issubset(mdata.mod):
            raise ValueError("checkpoint requires gene and guide modalities")
        if not mdata["gene"].obs_names.equals(mdata["guide"].obs_names):
            raise ValueError("gene and guide cell indexes differ")
        n_cells = int(mdata["gene"].n_obs)
        n_genes = int(mdata["gene"].n_vars)
        n_guides = int(mdata["guide"].n_vars)
    finally:
        mdata.file.close()
    manifest = {
        "status": "READY_FOR_INFERENCE",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_name": args.run_name,
        "pipeline_commit": args.commit,
        "mudata_filename": args.output.name,
        "mudata_sha256": sha256(args.output),
        "gtf_sha256": sha256(args.gtf),
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_guides": n_guides,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
