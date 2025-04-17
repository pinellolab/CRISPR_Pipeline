#!/usr/bin/env python

import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser(description='Process paired-end read files.')
    parser.add_argument('--dir', required=True, help='directory of fastq files')
    parser.add_argument('--reads1', required=True, help='list of read 1 files')
    parser.add_argument('--reads2', required=True, help='list of read 2 files')
    args = parser.parse_args()

    reads_1 = re.sub(r'\s', '', args.reads1)
    reads_2 = re.sub(r'\s', '', args.reads2)

    result = ' '.join(f'{os.path.join(args.dir, r1)} {os.path.join(args.dir, r2)}' for r1, r2 in zip(reads_1.split(';'), reads_2.split(';')))
    
    print(result)

if __name__ == '__main__':
    main()
