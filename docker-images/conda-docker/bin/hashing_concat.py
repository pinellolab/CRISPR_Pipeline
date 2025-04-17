#!/usr/bin/env python

import anndata as ad
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Process AnnData files and add batch information.')
    parser.add_argument('file_path_list', nargs='+', help='List of AnnData file paths')
    parser.add_argument('--output', type=str, required=True, help='Output file name')
    parser.add_argument('--temp_dir', type=str, default='./temp_processed', help='Directory for processed files')

    args = parser.parse_args()
    
    # Create temp dir if it doesn't exist
    os.makedirs(args.temp_dir, exist_ok=True)
    
    # Sort file paths
    sorted_file_paths = sorted(args.file_path_list, key=os.path.basename)
    print("Sorted file paths:", sorted_file_paths)
    
    var_index_name = None  # Store the var index name from the first file
    
    # Process and save each file
    processed_files = []
    for idx, file_path in enumerate(sorted_file_paths):
        print(f"Processing {file_path}")
        
        adata = ad.read_h5ad(file_path)
        
        if idx == 0 and adata.var_names.name is not None:
            var_index_name = adata.var_names.name
            print(f"Captured var index name: {var_index_name}")

        processed_file_path = os.path.join(args.temp_dir, f"processed_{idx}.h5ad")
        print(f"Saving processed file to {processed_file_path}")
        adata.write_h5ad(processed_file_path)
        processed_files.append(processed_file_path)
    
    # Use concat_on_disk with the processed file paths
    print(f"Concatenating files: {processed_files}")
    combined_adata = ad.experimental.concat_on_disk(
        processed_files, 
        join='outer', 
        out_file=args.output
    )
    
    if var_index_name:
        print(f"Restoring var index name: {var_index_name}")
        final_adata = ad.read_h5ad(args.output)
        final_adata.var_names.name = var_index_name
        final_adata.write_h5ad(args.output)
    
    print(f"Combined AnnData saved to {args.output}")

if __name__ == "__main__":
    main()