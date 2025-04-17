#!/usr/bin/env python

import argparse
import pandas as pd
import muon as mu

def select_nodes(mdata, num_nodes: int):
    results_df = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})
    results_df_sorted = results_df.sort_values(by='p_value')
    unique_intended_names = results_df_sorted['intended_target_name'].unique()
    selected = unique_intended_names[:num_nodes]
    selected_nodes = ' '.join(selected)
    
    return selected_nodes

def main():

    parser = argparse.ArgumentParser(description="Select nodes from MuData based on p_value ranking.")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData (.h5mu) file")
    parser.add_argument("--num_nodes", type=int, required=True, help="Number of central nodes to select")
    
    args = parser.parse_args()
    mdata = mu.read(args.mdata_path)
    selected_nodes = select_nodes(mdata, args.num_nodes)
    print(selected_nodes)

if __name__ == "__main__":
    main()
