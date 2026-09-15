#!/usr/bin/env python3
"""Report YES/NO if mudata.uns contains inference result keys."""
from __future__ import annotations

import argparse

from qc_mudata_io import result_keys

RESULT_KEYS = {
    "global_analysis_per_guide_results",
    "local_analysis_per_guide_results",
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

    # Only the key names decide this, so read the names and none of the tables:
    # loading them cost a 6.4 GB read to answer YES/NO.
    keys = set(result_keys(args.input))
    has_results = any(k in keys for k in RESULT_KEYS)
    print("YES" if has_results else "NO")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
