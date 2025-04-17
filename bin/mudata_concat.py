#!/usr/bin/env python
import argparse
import glob
import os
import mudata as md
import scipy.sparse as sp
import pandas as pd
import numpy as np
import copy
import gc

def concat_mudatas(input_pattern, output_file):
    """
    Concatenate multiple MuData files preserving all metadata, guide_assignment layers,
    and ensuring var DataFrames and uns dictionaries are properly maintained.
    """
    # Find all matching files
    files = sorted(glob.glob(input_pattern))
    if not files:
        print(f"No files found matching: {input_pattern}")
        return
    
    print(f"Found {len(files)} files to concatenate")
    
    print(f"Loading reference file: {os.path.basename(files[0])}")
    reference_mdata = md.read(files[0])
    
    # Create empty dictionaries to store metadata
    var_dataframes = {}
    uns_data = {}
    
    # Initialize from reference data
    for mod_name in reference_mdata.mod:
        var_dataframes[mod_name] = []
        uns_data[mod_name] = copy.deepcopy(reference_mdata[mod_name].uns)
    
    # Also keep global uns data
    global_uns = copy.deepcopy(reference_mdata.uns)
    
    # Check if there's only one file - handle specially
    if len(files) == 1:
        print("Only one file found, skipping concatenation but still processing metadata...")
        combined_mdata = reference_mdata
    else:
        mudatas = [reference_mdata] 
        
        for i, file in enumerate(files[1:], 1): 
            print(f"Loading {i+1}/{len(files)}: {os.path.basename(file)}")
            mdata = md.read(file)
            
            # Collect var DataFrames from each modality
            for mod_name in mdata.mod:
                if mod_name in var_dataframes:
                    # Make a copy of the var DataFrame with index as a regular column
                    mod_var = mdata[mod_name].var.copy()
                    mod_var['_var_name'] = mod_var.index
                    var_dataframes[mod_name].append(mod_var)
                    
                    # Merge any new uns data
                    if hasattr(mdata[mod_name], 'uns') and mdata[mod_name].uns:
                        for key, value in mdata[mod_name].uns.items():
                            if key not in uns_data[mod_name]:
                                uns_data[mod_name][key] = value
            
            # Also merge global uns data
            if hasattr(mdata, 'uns') and mdata.uns:
                for key, value in mdata.uns.items():
                    if key not in global_uns:
                        global_uns[key] = value
            
            mudatas.append(mdata)
        
        # Concatenate all MuData objects
        print("Concatenating all MuData objects...")
        combined_mdata = md.concat(mudatas)
        
        # Clean up individual MuData objects to free memory
        del mudatas
        gc.collect()
    
    # Apply preserved global uns data
    print("Restoring global uns data...")
    for key, value in global_uns.items():
        combined_mdata.uns[key] = value
    
    # Process each modality
    for mod_name in var_dataframes:
        print(f"Processing {mod_name} modality...")
        
        # Restore uns data for this modality
        print(f"  Restoring uns data for {mod_name}...")
        for key, value in uns_data[mod_name].items():
            combined_mdata[mod_name].uns[key] = value
        
        # Restore var data for this modality
        print(f"  Preserving var information for {mod_name}...")
        if var_dataframes[mod_name]:
            # Combine var DataFrames for this modality
            all_vars = pd.concat(var_dataframes[mod_name], axis=0)
            all_vars = all_vars.drop_duplicates(subset=['_var_name'])
            
            # Set index back to var names
            all_vars = all_vars.set_index('_var_name')
            current_vars = combined_mdata[mod_name].var
            matching_vars = all_vars.index.intersection(current_vars.index)
            
            if len(matching_vars) < len(current_vars.index):
                print(f"  Warning: Some vars in the combined object not found in collected data")
                print(f"  Found {len(matching_vars)} out of {len(current_vars.index)} vars")
                
                # Add missing indices to all_vars with empty data
                missing_vars = current_vars.index.difference(all_vars.index)
                empty_df = pd.DataFrame(index=missing_vars)
                empty_df.index.name = '_var_name'
                all_vars = pd.concat([all_vars, empty_df])
        
            all_vars = all_vars.loc[current_vars.index]
            
            for col in all_vars.columns:
                if col == '_var_name':
                    continue
                combined_mdata[mod_name].var[col] = all_vars[col]
    
    # Special handling for guide_assignment layer
    if 'guide' in combined_mdata.mod and 'guide_assignment' in reference_mdata['guide'].layers:
        print("Verifying guide_assignment layer...")
        
        # Ensure guide_assignment exists in the combined object
        if 'guide_assignment' not in combined_mdata['guide'].layers:
            print("Adding missing guide_assignment layer...")
            # We need to recreate it from scratch
            
            # Load guide_assignment matrices separately
            assignments = []
            for file in files:
                mdata = md.read(file)
                if 'guide' in mdata.mod and 'guide_assignment' in mdata['guide'].layers:
                    assignments.append(mdata['guide'].layers['guide_assignment'])
                del mdata
            
            # Stack the matrices if we found assignments
            if assignments:
                combined_matrix = sp.vstack(assignments)
                
                # Add the layer
                expected_shape = (combined_mdata.n_obs, combined_mdata['guide'].n_vars)
                if combined_matrix.shape == expected_shape:
                    combined_mdata['guide'].layers['guide_assignment'] = combined_matrix
                else:
                    print(f"Warning: guide_assignment shape mismatch: got {combined_matrix.shape}, expected {expected_shape}")
        else:
            # Verify the existing layer has the right shape
            current_shape = combined_mdata['guide'].layers['guide_assignment'].shape
            expected_shape = (combined_mdata.n_obs, combined_mdata['guide'].n_vars)
            
            if current_shape != expected_shape:
                print(f"Fixing guide_assignment layer shape: {current_shape} -> {expected_shape}")
                
                # Reload and stack matrices to get the correct shape
                assignments = []
                for file in files:
                    mdata = md.read(file)
                    if 'guide' in mdata.mod and 'guide_assignment' in mdata['guide'].layers:
                        assignments.append(mdata['guide'].layers['guide_assignment'])
                    del mdata
                
                if assignments:
                    combined_matrix = sp.vstack(assignments)
                    if combined_matrix.shape == expected_shape:
                        combined_mdata['guide'].layers['guide_assignment'] = combined_matrix
                    else:
                        print(f"Warning: Reconstructed guide_assignment shape still incorrect: {combined_matrix.shape}")
    
    # Save the combined result
    print(f"Saving combined MuData with {combined_mdata.n_obs} cells to {output_file}")
    combined_mdata.update()
    combined_mdata.write(output_file)
    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Concatenate MuData files")
    parser.add_argument("-i", "--input", dest="input", required=True, help="Input pattern (e.g., '*_output.h5mu')")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output file path")
    args = parser.parse_args()
    
    concat_mudatas(args.input, args.output)

if __name__ == "__main__":
    main()