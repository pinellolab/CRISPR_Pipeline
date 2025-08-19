#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
import numpy as np
import os


def system_call(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = result.stdout
    print(output)
    error = result.stderr
    return output, error

def process_reads(yaml_file):
    read_ids = []
    with open(yaml_file, 'r') as file:
        for line in file:
            if 'read_id:' in line:
                read_id = line.split('read_id:')[1].strip()
                read_ids.append(read_id)

    return ','.join(read_ids)


def get_info_from_yaml(modality, yaml_file, whitelist):
    reads = process_reads(yaml_file)
    print('reads out', reads)
    cmd = f"seqspec index -m {modality} -t kb -i {reads} {yaml_file}"
    print(cmd)
    output, error = system_call(cmd)

    # Print the output and error for debugging
    print("Output:", output)
    print("Error:", error)

    # add columns
    return pd.DataFrame([[modality, output.replace('\n', ''), whitelist, yaml_file]],
                        columns=['modality', 'representation', 'barcode_whitelist', 'seqspec_file'])


def main(modalities, yaml_file, output_file, whitelist):
    df_list = [get_info_from_yaml(m, yaml_file, whitelist) for m in modalities]
    pd.concat(df_list).to_csv(output_file, index=False, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process modalities and YAML file")
    parser.add_argument('--modalities', nargs='+', required=True, help="List of modalities")
    parser.add_argument('--yaml_file', required=True, help="Path to the YAML file")
    parser.add_argument('--output_file', required=True, help="Output file path")
    parser.add_argument('--whitelist', required=True, help="Path to the whitelist file")
    args = parser.parse_args()

    main(args.modalities, args.yaml_file, args.output_file, args.whitelist)
