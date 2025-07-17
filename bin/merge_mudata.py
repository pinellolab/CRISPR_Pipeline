#!/usr/bin/env python
import argparse
import muon as mu
import pandas as pd
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Combine cis and trans MuData files into a single MuData with separate test results")
    parser.add_argument("--cis_mudata", required=True, type=str,
                       help="Path to the cis MuData file (.h5mu)")
    parser.add_argument("--trans_mudata", required=True, type=str,
                       help="Path to the trans MuData file (.h5mu)")
    parser.add_argument("--output", required=True, type=str,
                       help="Output path for the combined MuData file (.h5mu)")
    parser.add_argument("--cis_results_key", default="test_results", type=str,
                       help="Key for test results in cis MuData (default: test_results)")
    parser.add_argument("--trans_results_key", default="test_results", type=str,
                       help="Key for test results in trans MuData (default: test_results)")

    args = parser.parse_args()

    # Validate input files exist
    if not os.path.exists(args.cis_mudata):
        print(f"Error: Cis MuData file not found: {args.cis_mudata}")
        sys.exit(1)

    if not os.path.exists(args.trans_mudata):
        print(f"Error: Trans MuData file not found: {args.trans_mudata}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        # Load MuData files
        print("Loading cis MuData file...")
        mudata_cis = mu.read(args.cis_mudata)

        print("Loading trans MuData file...")
        mudata_trans = mu.read(args.trans_mudata)

        # Check if the required keys exist
        if args.cis_results_key not in mudata_cis.uns:
            print(f"Error: Key '{args.cis_results_key}' not found in cis MuData")
            sys.exit(1)

        if args.trans_results_key not in mudata_trans.uns:
            print(f"Error: Key '{args.trans_results_key}' not found in trans MuData")
            sys.exit(1)

        # Extract cis results and add to trans MuData
        print("Extracting cis results...")
        cis_result = pd.DataFrame(mudata_cis.uns[args.cis_results_key])
        mudata_trans.uns['cis_test_results'] = cis_result

        # Rename trans results
        print("Renaming trans results...")
        mudata_trans.uns['trans_test_results'] = mudata_trans.uns.pop(args.trans_results_key)

        # Save the combined MuData
        print(f"Saving combined MuData to {args.output}...")
        mudata_trans.write(args.output)

        print("Successfully combined MuData files!")
        print(f"Output contains:")
        print(f"  - cis_test_results: {len(cis_result)} entries")
        print(f"  - trans_test_results: {len(mudata_trans.uns['trans_test_results'])} entries")

    except Exception as e:
        print(f"Error processing MuData files: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
