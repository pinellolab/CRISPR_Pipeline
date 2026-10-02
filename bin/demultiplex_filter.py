#!/usr/bin/env python

import anndata as ad
import pandas as pd
import argparse
import json

def filter_adata_demux(adata_path, demux_report_path, demux_config_path, filtered_output_path, unfiltered_output_path, demux_qc_path=None):
    # Load the AnnData object
    adata = ad.read_h5ad(adata_path)
    demux_report = pd.read_csv(demux_report_path, index_col=0)
    demux_config = pd.read_csv(demux_config_path, header=None)
    
    # Configure column names and data
    demux_config.columns = ['cluster_id', 'hto_type']
    demux_config['hto_type'] = demux_config['hto_type'].str.strip()
    adata.obs['cluster_id'] = demux_report['Cluster_id']
    adata.obs['gmm_demux_confidence'] = demux_report['Confidence']
    
    # Process multiplets
    demux_config['hto_type_split'] = demux_config['hto_type'].str.split('-').str.join(",")
    demux_config['hto_type_split'] = demux_config['hto_type_split'].apply(
        lambda x: 'multiplets' if ',' in x else x
    )
    
    # Merge demux config with adata.obs
    adata.obs = adata.obs.merge(
        demux_config,
        how='left',
        left_on='cluster_id',
        right_on='cluster_id'
    ).set_index(adata.obs.index)

    if demux_qc_path:
        with open(demux_qc_path) as handle:
            demux_qc = json.load(handle)
        # Older AnnData versions in the GMM-Demux image cannot serialize a
        # list of nested dictionaries in ``uns``.  Keep the full audit as a
        # canonical JSON string and expose the commonly queried values as
        # scalar fields.
        adata.uns['gmm_demux_qc_json'] = json.dumps(
            demux_qc, sort_keys=True, separators=(',', ':')
        )
        adata.uns['gmm_demux_selected_seed'] = int(demux_qc['selected_seed'])
        adata.uns['gmm_demux_attempt_count'] = int(len(demux_qc['attempts']))
    
    # Filter data
    adata_filtered = adata[
        (adata.obs['hto_type'] != 'negative') & 
        (adata.obs['hto_type_split'] != 'multiplets')
    ].copy()
    
    # Save outputs
    adata.write(unfiltered_output_path)
    adata_filtered.write(filtered_output_path)

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Update AnnData object with demultiplexing report."
    )
    parser.add_argument(
        '--demux_qc',
        default=None,
        help="JSON report from deterministic GMM-Demux seed selection."
    )
    parser.add_argument(
        '--adata', 
        type=str, 
        required=True, 
        help="Path to the input h5ad file."
    )
    parser.add_argument(
        '--demux_report', 
        type=str, 
        required=True, 
        help="Path to the demultiplexing report CSV file."
    )
    parser.add_argument(
        '--demux_config', 
        type=str, 
        required=True, 
        help="Path to the demultiplexing config CSV file."
    )
    parser.add_argument(
        '--filtered_output', 
        type=str, 
        required=True, 
        help="Path to save the filtered h5ad file."
    )
    parser.add_argument(
        '--unfiltered_output', 
        type=str, 
        required=True, 
        help="Path to save the unfiltered h5ad file."
    )
    
    args = parser.parse_args()
    filter_adata_demux(
        args.adata, 
        args.demux_report, 
        args.demux_config, 
        args.filtered_output,
        args.unfiltered_output,
        args.demux_qc
    )

if __name__ == "__main__":
    main()
