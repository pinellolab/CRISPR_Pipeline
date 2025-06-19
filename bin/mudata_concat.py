#!/usr/bin/env python
import argparse
import os
import mudata as md

def concat_mudatas(input_files, output_file):
    """
    Concatenate multiple MuData files. If only one file is provided, it's copied to the output.
    """
    files = sorted(input_files, key=lambda x: os.path.basename(x))
    if not files:
        print(f"No files found: {input_files}")
        return

    print(f"Found {len(files)} files to concatenate")

    # Handle single file case
    if len(files) == 1:
        print(f"Only one file found. Copying {files[0]} to {output_file}")
        single_mdata = md.read(files[0])
        print(f"Saving MuData with {single_mdata.n_obs} cells to {output_file}")
        single_mdata.write(output_file)
        return

    # Handle multiple files case
    print("Concatenating all MuData objects...")
    combined_mdata = md.concat([md.read(ff) for ff in files], merge='first', uns_merge='first', join='outer')

    print(f"Saving combined MuData with {combined_mdata.n_obs} cells to {output_file}")
    combined_mdata.write(output_file)

    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Concatenate MuData files")
    parser.add_argument("-i", "--input", dest="input", nargs="+", required=True, help="Input mudata files")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output file path")
    args = parser.parse_args()

    concat_mudatas(args.input, args.output)

if __name__ == "__main__":
    main()
