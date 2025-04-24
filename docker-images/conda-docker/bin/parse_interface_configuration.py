#!/usr/bin/env python
import argparse
import pandas as pd
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Generate configuration files for sequencing data.')
    
    # Required arguments
    parser.add_argument('--config_table', required=True, help='Path to the configuration table CSV file')
    parser.add_argument('--output_dir', required=True, help='Directory to save the output configuration files')

    # Optional arguments
    parser.add_argument('--mapping_tabs', nargs='+', default=['scRNA', 'Guides', 'Hash'], 
                        help='List of mapping tabs to process (default: scRNA, Guides, Hash)')
    
    return parser.parse_args()

def format_fastq_files(grouped_df, tab_name):
    grouped = grouped_df.groupby('batch_name').apply(
        lambda x: ' '.join([f"{r1} {r2}" for r1, r2 in zip(x['read1'], x['read2'])])
    )
    
    # Format according to tab name
    if tab_name == 'scRNA':
        config_key = 'fastq_files_rna'
    elif tab_name == 'Guides':
        config_key = 'fastq_files_guide'
    elif tab_name == 'Hash':
        config_key = 'fastq_files_hashing'
    else:
        config_key = 'fastq_files'
    
    return f'{config_key} = [\n    ' + ',\n    '.join(f'"{seq}"' for seq in grouped) + '\n    ]'

def format_test_fastq_files(filtered_df, tab_name):
    filtered_df = filtered_df.dropna(subset=['read1', 'read2'])
    read1_list = filtered_df['read1'].apply(str.strip).tolist()
    read2_list = filtered_df['read2'].apply(str.strip).tolist()

    if tab_name == 'Guides':
        r1_key = 'test_guide_fastq_r1'
        r2_key = 'test_guide_fastq_r2'
    elif tab_name == 'Hash':
        r1_key = 'test_hashing_fastq_r1'
        r2_key = 'test_hashing_fastq_r2'
    else:
        return ""

    return f'{r1_key} = {read1_list}\n    {r2_key} = {read2_list}'

def generate_config_files(config_table_path, output_dir, mapping_tabs):
    config_table = pd.read_csv(config_table_path)

    ### Write sequence parameters
    config_sequences = {}
    dataset_hashing_status = 'false'
    test_fastq_configs = []

    for tab_name in mapping_tabs:
        filtered_df = config_table[config_table['tab_name'] == tab_name]
        if not filtered_df.empty:
            config_sequences[tab_name] = format_fastq_files(filtered_df, tab_name)
            if tab_name in ['Guides', 'Hash']:
                test_fastq_configs.append(format_test_fastq_files(filtered_df, tab_name))
            if tab_name == 'Hash':
                dataset_hashing_status = 'true'

    ### Write non-sequence parameters
    config_table_nonseq = config_table[config_table['variable'] != 'sequence']
    config_table_nonseq2 = config_table_nonseq[['variable', 'variable_value']]

    config_nonseq = []
    distance_from_center_present = False

    for _, row in config_table_nonseq2.iterrows():
        var_name = row['variable']
        var_value = row['variable_value']

        if var_name == 'distance_from_center':
            distance_from_center_present = True
        
        # Check if the value is numeric 
        try:
            if '.' in var_value:
                var_value = float(var_value)
            else:
                var_value = int(var_value)
        except ValueError:
            var_value = f"'{var_value}'"
        
        config_nonseq.append(f"{var_name} = {var_value}")
    
    # Add distance_from_center if not present
    if not distance_from_center_present:
        config_nonseq.append("distance_from_center = 1000000")

    config_nonseqs = '\n    '.join(config_nonseq)

    ### Write covariate list
    df_non_nan = config_table[config_table['batch_name'].notna()]
    split_batches = df_non_nan['batch_name'].apply(lambda x: x.split(', '))

    max_splits = split_batches.apply(len).max()
    column_names = ['batch'] + [f'cov{i}' for i in range(1, max_splits)]

    batch = split_batches.apply(lambda x: x[0])
    covariates = pd.DataFrame(split_batches.tolist(), columns=column_names)
    covariates = covariates.drop_duplicates()

    covariate_dict = covariates.drop(columns=['batch']).to_dict(orient='list')

    covariate_list = f"params.covariate_list = [\n    batch: {batch.unique().tolist()},\n"
    batch_only = f"batch={batch.unique().tolist()}"
    for cov_name, cov_values in covariate_dict.items():
        covariate_list += f"    {cov_name}: {cov_values},\n"

    covariate_list = covariate_list.rstrip(',\n') + '\n]'

    # Combine all config sections inside params
    config_sequences_combined = ''.join([f'{config}\n    ' for config in config_sequences.values()])
    test_fastq_configs_combined = '\n\n    '.join(test_fastq_configs)

    final_config_output = (
        f"params {{\n"
        f"    DATASET_HASHING = '{dataset_hashing_status}'\n\n"
        f"    {config_nonseqs}\n\n"
        f"    {config_sequences_combined}\n\n"
        f"    {test_fastq_configs_combined}\n"
        f"    {batch_only}\n"
        f"}}\n\n"
        f"{covariate_list}"
    )

    # Write final config
    final_config_file_path = os.path.join(output_dir, 'pipeline_input.config')
    with open(final_config_file_path, 'w') as file:
        file.write(final_config_output)

    print(f"Configuration files saved as pipeline_input.config. It can be found in {output_dir}")

def main():
    args = parse_args()

    # Generate config files using provided arguments
    generate_config_files(
        config_table_path=args.config_table,
        output_dir=args.output_dir,
        mapping_tabs=args.mapping_tabs
    )

if __name__ == '__main__':
    main()
