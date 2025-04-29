#!/usr/bin/env python
import argparse
import glob
import os
import mudata as md
import scipy.sparse as sp
import pandas as pd
import numpy as np
import copy
import gc

def concat_mudatas(input_files, output_file):
    """
    Concatenate multiple MuData files preserving all metadata, guide_assignment layers,
    and ensuring var DataFrames and uns dictionaries are properly maintained.
    """
    # Sort files so that batches are grouped together (ignoring subfolders if present)
    files = sorted(input_files, key = lambda x: os.path.basename(x))
    if not files:
        print(f"No files found: {input_files}")
        return
    
    print(f"Found {len(files)} files to concatenate")
    

    # Concatenate all MuData objects
    print("Concatenating all MuData objects...")
    combined_mdata = md.concat([md.read(ff) for ff in files], 
        merge = 'first', uns_merge = 'first', join = 'outer')

    # Save the combined result
    print(f"Saving combined MuData with {combined_mdata.n_obs} cells to {output_file}")
    combined_mdata.update()
    combined_mdata.write(output_file)
    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Concatenate MuData files")
    parser.add_argument("-i", "--input", dest="input", nargs="+", required=True, help="Input mudata files to concatenate")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output file path")
    args = parser.parse_args()
    
    concat_mudatas(args.input, args.output)

if __name__ == "__main__":
    main()