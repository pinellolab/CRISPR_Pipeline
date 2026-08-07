#!/usr/bin/env python
"""Convert an inference result table between supported serializations."""

import argparse
from pathlib import Path

from result_table_io import read_result_table, write_result_table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--parquet-compression",
        default="zstd",
        choices=("zstd", "snappy", "gzip", "none"),
    )
    args = parser.parse_args()

    if not args.input.is_file():
        parser.error(f"input result table does not exist: {args.input}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    frame = read_result_table(args.input)
    compression = (
        None if args.parquet_compression == "none" else args.parquet_compression
    )
    write_result_table(frame, args.output, parquet_compression=compression)
    print(
        f"Converted {len(frame)} rows: {args.input} -> {args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
