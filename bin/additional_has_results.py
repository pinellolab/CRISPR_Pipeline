#!/usr/bin/env python3
"""Report YES/NO if mudata.uns contains inference result keys."""
from __future__ import annotations

import argparse

import mudata

RESULT_KEYS = {
    "trans_per_guide_results",
    "per_guide_results",
    "cis_per_guide_results",
    "trans_test_results",
    "test_results",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    args = parser.parse_args()

    m = mudata.read_h5mu(args.input)
    keys = set(m.uns.keys())
    has_results = any(k in keys for k in RESULT_KEYS)
    print("YES" if has_results else "NO")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
