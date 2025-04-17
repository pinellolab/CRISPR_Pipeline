#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf
import muon as mu

def main(limit, mudata_path, input_gtf):
    mudata = mu.read(mudata_path)
    guide_data = mudata.mod['guide'].var
    
    # Load GTF
    df_gtf_refseq = read_gtf(input_gtf).to_pandas()
    df_gtf_unique_gene_name = df_gtf_refseq.drop_duplicates('gene_name')
    df_gtf_unique_gene_name.set_index('gene_name', inplace=True)
    
    final_candidates_list = []
    for g_name, target_name, guide_chr, guide_start in guide_data[['guide_id', 'intended_target_name', 'guide_chr', 'guide_start']].values:
        if limit is None:
            # Use all genes from mudata for unlimited mode
            gene_df = pd.DataFrame({'guide_id': g_name,
                                    'gene_id': list(mudata['gene'].var_names),
                                    'intended_target_name': target_name,
                                    'pair_type': 'discovery'})
            final_candidates_list.append(gene_df)
        else:
            # Use distance-based filtering
            temp_gtf = df_gtf_unique_gene_name.query(f'seqname == "{guide_chr}"').copy()
            temp_gtf['distance_from_guide'] = np.abs(temp_gtf['start'] - guide_start)
            final_candidates = temp_gtf.query(f'distance_from_guide <= {limit}').copy()
            final_candidates['guide_id'] = g_name
            final_candidates['intended_target_name'] = target_name
            final_candidates['pair_type'] = 'discovery'
            final_candidates.reset_index(inplace=True)
            final_candidates_list.append(final_candidates[['guide_id', 'gene_id', 'intended_target_name', 'pair_type']].copy())
    
    final_candidates_df = pd.concat(final_candidates_list, ignore_index=True)
    final_candidates_df['gene_id'] = final_candidates_df['gene_id'].str.split('.').str[0]
    final_candidates_df = final_candidates_df.rename(columns={'gene_id': 'gene_name'})
    final_candidates_df.to_csv('pairs_to_test.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GTF and guide data.")
    parser.add_argument('--limit', type=int, default=1000000, help='Limit for distance from guide. Use -1 for no limit.')
    parser.add_argument('mudata', type=str, help='Path to the input MuData file')
    parser.add_argument('input_gtf', type=str, help='Path to the input GTF file')
    
    args = parser.parse_args()
    limit = None if args.limit == -1 else args.limit
    main(limit, args.mudata, args.input_gtf)