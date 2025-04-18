#!/usr/bin/env python

import pandas as pd
import argparse

def process_table(guide_table):
    guide_metadata = pd.read_csv(guide_table, sep="\t")
    guide_metadata[['spacer', 'guide_id']].to_csv('guide_features.txt',
                                                           sep='\t', header=None, index=None)

def main():
    parser = argparse.ArgumentParser(description="Process guide metadata.")
    parser.add_argument('--guide_table', type=str, required=True, help="Path to the input Guide file.")
    args = parser.parse_args()

    process_table(args.guide_table)

if __name__ == "__main__":
    main()
