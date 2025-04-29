#!/usr/bin/env python

import anndata as ad
import pandas as pd
import scipy.io
import gzip
import os
import argparse

def prepare_demux(input_file, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read in adata
    adata = ad.read_h5ad(input_file)

    # Parse the HTO tags
    var_names_list = adata.var_names.to_list()
    var_names_string = ",".join(var_names_list)

    # Exporting adata.obs as barcode.tsv.gz
    obs = pd.DataFrame(adata.obs_names.to_list())
    obs.to_csv(os.path.join(output_dir, 'barcodes.tsv'), index=False, header=False) 
    with open(os.path.join(output_dir, 'barcodes.tsv'), 'rb') as f_in, gzip.open(os.path.join(output_dir, 'barcodes.tsv.gz'), 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove(os.path.join(output_dir, 'barcodes.tsv'))  

    # Exporting adata.var as features.tsv.gz
    var = adata.var
    var.to_csv(os.path.join(output_dir, 'features.tsv'), sep='\t', index=True, header=False)  
    with open(os.path.join(output_dir, 'features.tsv'), 'rb') as f_in, gzip.open(os.path.join(output_dir, 'features.tsv.gz'), 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove(os.path.join(output_dir, 'features.tsv')) 

    # Exporting adata.X as matrix.mtx.gz
    scipy.io.mmwrite(os.path.join(output_dir, 'matrix.mtx'), adata.X.T) 
    with open(os.path.join(output_dir, 'matrix.mtx'), 'rb') as f_in, gzip.open(os.path.join(output_dir, 'matrix.mtx.gz'), 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove(os.path.join(output_dir, 'matrix.mtx')) 

    return var_names_string

def main():
    # Set up argparse 
    parser = argparse.ArgumentParser(description="Prepare script to export data from an AnnData object.")
    parser.add_argument('--adata', type=str, required=True, help="Path to the input h5ad file.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Path to the output directory where the files will be saved.")
    args = parser.parse_args()

    var_names_string = prepare_demux(args.adata, args.output)
    print(var_names_string)

if __name__ == "__main__":
    main()
