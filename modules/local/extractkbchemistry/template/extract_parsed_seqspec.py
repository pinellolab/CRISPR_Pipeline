#!/usr/bin/env python

import argparse
import os

import pandas as pd


def process_files(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Cannot find the parsed_seqspec file: {file_path}")

    df = pd.read_csv(file_path, sep='\t')
    return ''.join(df['representation'].astype(str))


def main():
    parser = argparse.ArgumentParser(description="Process parsed_seqspec files and extract seqspec info")
    parser.add_argument('--file', required=True, help="Path to the parsed_seqspec file.")
    args = parser.parse_args()

    print(process_files(args.file))


if __name__ == "__main__":
    main()