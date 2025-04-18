#!/usr/bin/env python

import argparse
import pandas as pd
import muon as mu

def main(guide_inference, mudata_path):
    # read in files
    print(f"Reading guide inference file from {guide_inference}...")
    guide_inference = pd.read_csv(guide_inference)
    print(f"Reading mudata file from {mudata_path}...")
    mudata = mu.read_h5mu(mudata_path)

    # Check if mudata.mod contains 'gene' and 'guide'
    if 'gene' not in mudata.mod or 'guide' not in mudata.mod:
        raise KeyError("Mudata file is missing 'gene' or 'guide' modality.")
    
    # Debugging: Print basic information about the input data
    print(f"Tested pairs dataset contains {len(guide_inference)} rows.")
    print(f"Mudata 'gene' modality contains {mudata.mod['gene'].shape[1]} genes.")
    print(f"Mudata 'guide' modality contains {mudata.mod['guide'].shape[1]} guides.")
    
    gene_var = mudata.mod['gene'].var.index if mudata.mod['gene'].var.index is not None else []
    guide_var = mudata.mod['guide'].var['guide_id'] if mudata.mod['guide'] is not None else []

    include1 = set(guide_inference['gene_name']).intersection(set(gene_var))
    include2 = set(guide_inference['guide_id']).intersection(set(guide_var))
    
    print(f"Number of genes in common: {len(include1)}")
    print(f"Number of guides in common: {len(include2)}")
    
    subset = guide_inference[guide_inference['gene_name'].isin(include1) & guide_inference['guide_id'].isin(include2)]
    
    # Ensure that subset is not empty
    if subset.empty:
        raise ValueError("The subset of guide_inference is empty after filtering. Please check your input data.")
    
    print(f"Subset contains {len(subset)} rows after filtering.")

    mudata.uns['pairs_to_test'] = subset.to_dict(orient='list')

    # rename keys in uns
    key_mapping = {
        'gene_name': 'gene_id',
        'guide_id': 'guide_id',
        'intended_target_name': 'intended_target_name',
        'pair_type': 'pair_type'
    }

    # Ensure pairs_to_test is not None before trying to access items
    if mudata.uns.get('pairs_to_test') is None:
        raise ValueError("'pairs_to_test' in mudata.uns is None, something went wrong when processing the subset.")

    # Update the keys in pairs_to_test
    print("Renaming keys in pairs_to_test...")
    mudata.uns['pairs_to_test'] = {key_mapping.get(k, k): v for k, v in mudata.uns['pairs_to_test'].items()}
    
    # save the mudata
    output_file = "mudata_inference_input.h5mu"
    print(f"Saving processed mudata to {output_file}...")
    mudata.write(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process mudata and add guide_inference to uns.")
    parser.add_argument('guide_inference', type=str, help='Path to the input inference file')
    parser.add_argument('mudata_path', type=str, help='Path to the input MuData file')

    args = parser.parse_args()
    main(args.guide_inference, args.mudata_path)
