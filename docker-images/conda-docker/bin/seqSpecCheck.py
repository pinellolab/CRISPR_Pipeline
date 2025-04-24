#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import gzip
from collections import Counter
import numpy as np  
import os

def readFastq(filename, max_reads):
    """Reads a FASTQ file and returns sequences, up to a maximum of max_reads."""
    sequences = []
    with gzip.open(filename, 'rt') as fh:
        while len(sequences) < max_reads:
            header = fh.readline()  # Read the header line
            if not header:
                break  # End of file
            seq = fh.readline().rstrip()  # Read base sequence
            fh.readline()  # Skip the '+' line
            fh.readline()  # Skip the quality line
            sequences.append(seq)
    return sequences

def fastq_sequence_plot(seqs, file_name, ax):
    """Plot nucleotide frequencies for a list of sequences on the given axis."""
    df_read = pd.DataFrame([list(n) for n in seqs])
    freq = []
    for n in ['A', 'C', 'G', 'T']:
        freq.append((df_read == n).sum() / len(df_read))
    freq_df = pd.DataFrame(freq, index=['A', 'C', 'G', 'T']).T
    freq_df.plot(color=['red', 'green', 'blue', 'orange'], ax=ax)
    ax.set_title(f'Nucleotide Frequency Plot - {file_name}')
    ax.set_xlabel('Position in Sequence')
    ax.set_ylabel('Frequency')

def find_sequence_positions(fastq_file, guide_metadata, max_reads=100000, sep='\t'):
    # Read FASTQ sequences
    seqs = readFastq(fastq_file, max_reads)
    
    # Load metadata and check for required column
    metadata_df = pd.read_csv(guide_metadata, sep=sep)
    guide_column = 'spacer' if 'spacer' in metadata_df else 'sequence'
    if guide_column not in metadata_df:
        raise ValueError("Guide metadata must contain 'spacer' or 'sequence' column.")
    
    positions = [
        s.index(n) for s in seqs for n in metadata_df[guide_column].values if n in s
    ]
    
    pos_table = (
            pd.DataFrame.from_dict(Counter(positions), orient='index', columns=['Count'])
            .rename_axis('Position Index')
            .sort_values(by='Count', ascending=False)
        )
    return pos_table

def main():
    parser = argparse.ArgumentParser(description="FASTQ sequence analysis tool")
    parser.add_argument('--read1', nargs='+', help="Paths to multiple gzipped FASTQ files for read1")
    parser.add_argument('--read2', nargs='+', help="Paths to multiple gzipped FASTQ files for read2")
    parser.add_argument('--metadata', help="Path to the metadata file")
    parser.add_argument('--max_reads', type=int, default=100000, help="Maximum number of reads to process")
    parser.add_argument('--plot', action='store_true', help="Plot the nucleotide frequencies")
    
    args = parser.parse_args()

    output_dir = 'seqSpec_plots'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Fix: Correct reference to args.read1 and args.read2
    read1 = args.read1
    read2 = args.read2
    guide_metadata = args.metadata
    max_reads = args.max_reads

    assert len(read1) == len(read2), "The number of read1 and read2 files must be the same."

    num_files = len(read1)

    # Create subplots with 2 columns: one for read1 and one for read2
    fig, axs = plt.subplots(nrows=num_files, ncols=2, figsize=(12, num_files * 3))

    # Ensure axs is always a 2D array
    if num_files == 1:
        axs = np.array([axs])  # Convert to 2D array for consistent indexing
    elif axs.ndim == 1:
        axs = axs.reshape(-1, 2) 

    all_r2_tables = []
    for i, (r1, r2) in enumerate(zip(read1, read2)):
        read1_name = os.path.basename(r1).replace('.fastq.gz', '')
        read2_name = os.path.basename(r2).replace('.fastq.gz', '')

        # Process read1
        sequences_r1 = readFastq(r1, max_reads)
        # Plot read1 sequences on axs[i, 0]
        fastq_sequence_plot(sequences_r1, read1_name, axs[i, 0])

        # Process read2
        sequences_r2 = readFastq(r2, max_reads)
        # Plot read2 sequences on axs[i, 1]
        fastq_sequence_plot(sequences_r2, read2_name, axs[i, 1])

        r2_table = find_sequence_positions(r2, guide_metadata, max_reads)
        read_name_row = pd.DataFrame({'Position Index': ['read name'], 'Count': [read2_name]})
        r2_table_with_name = pd.concat([read_name_row, r2_table.reset_index()], ignore_index=True)
        all_r2_tables.append(r2_table_with_name)

        # Save the final grid plot
    if args.plot:
        plt.tight_layout()
        plt.savefig("seqSpec_plots/seqSpec_check_plots.png", dpi=300)
        plt.close()

    # Concatenate all tables and save as a single CSV
    if all_r2_tables:
        final_r2_table = pd.concat(all_r2_tables, ignore_index=True)
        final_r2_table.to_csv("position_table.csv", index=False)

if __name__ == "__main__":
    main()
