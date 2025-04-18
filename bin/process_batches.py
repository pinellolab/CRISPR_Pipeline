#!/usr/bin/env python

import argparse
import os

def process_batches(fastq, dir):
    processed_batches = []
    batches = fastq[0].split(';')
    for batch in batches:
        processed_batch = batch.replace('|', ' ').replace(',', ' ')
        processed_batch = ' '.join(processed_batch.split())
        processed_batch = ' '.join(os.path.join(dir, filename) for filename in processed_batch.split())
  
        processed_batches.append(processed_batch)
    return processed_batches

def main():
    parser = argparse.ArgumentParser(description="Process fastq batches and prepend directory.")
    parser.add_argument('--fastq', type=str, nargs='+', required=True, help="List of fastq file batches.")
    parser.add_argument('--dir', type=str, required=True, help="Directory to prepend to each fastq file.")

    args = parser.parse_args()
    processed_batches = process_batches(args.fastq, args.dir)
    for i, batch in enumerate(processed_batches):
        print(f"{i} {batch}")

if __name__ == "__main__":
    main()
