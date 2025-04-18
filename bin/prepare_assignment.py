#!/usr/bin/env python
import argparse
import os
import mudata as md

def main():
    parser = argparse.ArgumentParser(description="Split MuData by batch")
    parser.add_argument("--input", help="Input MuData file (.h5mu)")
    parser.add_argument("-o", "--output-dir", default="./", help="Output directory")
    parser.add_argument("-b", "--batch-key", default="batch", help="Batch key in .obs")
    args = parser.parse_args()
    
    mdata = md.read(args.input)
    os.makedirs(args.output_dir, exist_ok=True)
    batches = mdata.obs[args.batch_key].unique()
    
    for batch in batches:
        batch_str = str(batch).replace(" ", "_")
        batch_data = mdata[mdata.obs[args.batch_key] == batch, :]
        
        batch_data.update()
        output_file = os.path.join(args.output_dir, f"{batch_str}_mudata.h5mu")
        batch_data.write(output_file)

if __name__ == "__main__":
    main()