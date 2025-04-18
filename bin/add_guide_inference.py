#!/usr/bin/env python
import muon as mu
import pandas as pd
import argparse

def add_guide_inference(test_results_csv, mudata_path):
    # Load test results and MuData
    test_results = pd.read_csv(test_results_csv)
    mudata = mu.read_h5mu(mudata_path)

    # RENAME 
    test_results.rename(columns={
        'log2_fc': 'sceptre_log2_fc',
        'p_value': 'sceptre_p_value'
    }, inplace=True)
    
    # Update MuData with test results
    mudata.uns["test_results"] = test_results.to_dict(orient='list')
    mudata.write("inference_mudata.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process test results and update MuData object.')
    parser.add_argument('--test_results_csv', required=True, help='Path to the test results CSV file.')
    parser.add_argument('--mudata', required=True, help='Path to the input h5mu file.')

    args = parser.parse_args()
    add_guide_inference(args.test_results_csv, args.mudata)
