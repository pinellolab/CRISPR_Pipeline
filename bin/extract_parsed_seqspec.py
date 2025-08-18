#!/usr/bin/env python

import pandas as pd
import argparse
import os

def process_files(file_path):
    if not os.path.isfile(file_path):
        print("Cannot find the parsed_seqspec file")
        return None

    df = pd.read_csv(file_path, sep='\t')
    return ''.join(df['representation'].astype(str))

def main():
    parser = argparse.ArgumentParser(description="Process parsed_seqspec files and extract seqspec info")
    parser.add_argument('--file', required=True, help="Path to the parsed_seqspec file.")

    args = parser.parse_args()

    result = process_files(args.file)
    if result:
        print(result)

if __name__ == "__main__":
    main()
