#!/usr/bin/env python

import os
import re
import shutil 
import argparse

def process_directory(base_dir, output_dir, pattern, file_name, prefix=''):
    match = pattern.match(file_name)
    if match:
        extracted_info = match.group(1)
        item_path = os.path.join(base_dir, file_name)
        inspect_json_path = os.path.join(item_path, "inspect.json")
        run_info_json_path = os.path.join(item_path, "run_info.json")
        if os.path.exists(inspect_json_path) and os.path.exists(run_info_json_path):
            new_inspect_json_path = os.path.join(output_dir, f"{prefix}-{extracted_info}-inspect.json")
            new_run_info_json_path = os.path.join(output_dir, f"{prefix}-{extracted_info}-run_info.json")
            shutil.copy2(inspect_json_path, new_inspect_json_path)
            shutil.copy2(run_info_json_path, new_run_info_json_path)
            print(f"Copied: {inspect_json_path} to {new_inspect_json_path}")
            print(f"Copied: {run_info_json_path} to {new_run_info_json_path}")
        else:
            print(f"JSON files not found in {item_path}")
    else:
        print(f"No match for: {file_name}")

def main(base_dir, output_dir):
    transcript_pattern = re.compile(r'(.+)_ks_transcripts_out')
    guide_pattern = re.compile(r'(.+)_ks_guide_out')

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            process_directory(base_dir, output_dir, transcript_pattern, item, prefix='trans')
            process_directory(base_dir, output_dir, guide_pattern, item, prefix='guide')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process directories and rename files based on patterns.")
    parser.add_argument(
        "base_dir", 
        type=str, 
        nargs='?',  
        default=os.getcwd(),  # Default to the current working directory
        help="The base directory containing the directories to process. Defaults to the current directory."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        nargs='?', 
        default=os.path.join(os.getcwd(), "json_dir"),  
        help="The directory where the renamed files will be saved. Defaults to a new 'output' directory in the current directory."
    )
    
    args = parser.parse_args()
    main(args.base_dir, args.output_dir)

