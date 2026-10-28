#!/usr/bin/env python3
"""Apply the configured action to a guide-mapping QC report."""

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--action", choices=("error", "warn"), default="error")
    args = parser.parse_args()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    if report.get("status") == "PASS":
        print("Guide mapping QC passed.")
        return 0
    message = "Guide mapping QC failed: " + " | ".join(report.get("failures", []))
    if args.action == "warn":
        print("WARNING: " + message, file=sys.stderr)
        return 0
    print("ERROR: " + message, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
