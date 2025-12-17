#!/usr/bin/env python
import pandas as pd
from Bio.Seq import Seq
import crispr_ambiguous_mapping
from importlib.metadata import version
import argparse
import os
import subprocess # Use subprocess for shell commands
import anndata
import numpy as np


import numpy as np

def extract_fasta_sequences(string_value):
  string_decripted = np.array(string_value.split(' ')).reshape(-1,2)
  R1, R2 = string_decripted[:,0].flatten().tolist(),  string_decripted[:,1].flatten().tolist()
  R1 = ' '.join(R1)
  R2 = ' '.join(R2)
  return R1, R2




def parse_seqspec_string(chemistry: str):
    """
    Strict parser for 'a:b:c,d:e:f,g:h:i' describing bc:umi:seq (in that order).
    Returns: protospacer_start, protospacer_len, umi_start, umi_len, barcode_start, barcode_len
    """
    s = chemistry.strip(':')
    print (s, 'strip')
    parts = [p.strip() for p in s.split(':')]
    print (parts)
    if len(parts) != 3:
        raise ValueError("Expected 3 triplets separated by commas: bc,umi,seq")

    def parse_triplet(tri: str):
        fields = [x.strip() for x in tri.split(',')]
        if len(fields) != 3:
            raise ValueError(f"Invalid triplet '{tri}'. Use 'file:start:stop'.")
        f, start, stop = (int(fields[0]), int(fields[1]), int(fields[2]))
        return f, start, stop

    (bc_f, bc_start, bc_stop)   = parse_triplet(parts[0])  # barcode
    (umi_f, umi_start, umi_stop)= parse_triplet(parts[1])  # UMI
    (seq_f, seq_start, seq_stop)= parse_triplet(parts[2])  # protospacer/seq

    def to_len(start, stop):
        return None if stop == 0 else (stop - start)

    barcode_len     = to_len(bc_start, bc_stop)
    umi_len         = to_len(umi_start, umi_stop)
    protospacer_len = to_len(seq_start, seq_stop)

    return (seq_start, protospacer_len, umi_start, umi_len, bc_start, barcode_len)

# Example
# chemistry = "0:0:16,0:16:28,1:63:83"
# protospacer_start, protospacer_len, umi_start, umi_len, barcode_start, barcode_len = parse_seqspec_string(chemistry)


def series_to_anndata(count_series: pd.Series) -> anndata.AnnData:
    """
    Convert a multi-indexed Series (CellBarcode, protospacer) to an AnnData object.

    Parameters
    ----------
    count_series : pd.Series
        Series with MultiIndex (CellBarcode, protospacer) and numeric values.

    Returns
    -------
    adata : anndata.AnnData
        AnnData object with cells as obs, protospacers as var, and counts as X.
    """
    # Convert Series to DataFrame with cells as rows, protospacers as columns
    df = count_series.unstack(level='protospacer').fillna(0)

    # Create AnnData
    adata = anndata.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    adata.var['guide_id'] = adata.var.index
    return adata

def run_crispr_mapping(args):
    """
    Runs the CRISPR guide mapping workflow using the provided arguments.
    """
    print(f"🔬 Starting CRISPR Mapping (crispr_ambiguous_mapping v{version('crispr_ambiguous_mapping')})")
    print("-" * 40)

    fastq_r1_fn, fastq_r2_fn = extract_fasta_sequences(args.fastq)

    # --- 1. Load Input Files ---
    
    print (args.guide_set_fn, ' this is the name')
    guide_set_df = pd.read_csv(str(args.guide_set_fn), sep='\t')
    print ('after pd.read_csv')
    barcode_inclusion_df = pd.read_csv(args.barcode_inclusion_list_fn, header=None)
    barcode_inclusion_df.columns = ["barcode"]

    
    guide_whitelist_input_df = pd.DataFrame(guide_set_df["spacer"])
    guide_whitelist_input_df.columns = ["protospacer"]
    print(f"Loaded {len(guide_whitelist_input_df)} unique protospacers from guide set.")
    #print ('forcing downsample REMOVE IT')
    #args.downsample_reads = 500_000
    # --- 2. Downsample FASTQ Files ---
    if args.downsample_reads > 0:
        print(f"Downsampling to {args.downsample_reads} read pairs...")
        fastq_r1_downsampled_fn = os.path.basename(fastq_r1_fn).replace(".fastq.gz", f".downsampled_{args.downsample_reads}.fastq")
        fastq_r2_downsampled_fn = os.path.basename(fastq_r2_fn).replace(".fastq.gz", f".downsampled_{args.downsample_reads}.fastq")

        read_count = args.downsample_reads * 4 # 4 lines per read entry (name, seq, +, qual)
        
        # Downsample R1
        subprocess.run(f"zcat {fastq_r1_fn} | head -n {read_count} > {fastq_r1_downsampled_fn}", shell=True, check=True, executable='/bin/bash')
        # Downsample R2
        subprocess.run(f"zcat {fastq_r2_fn} | head -n {read_count} > {fastq_r2_downsampled_fn}", shell=True, check=True, executable='/bin/bash')
        
        r1_input = [fastq_r1_downsampled_fn]
        r2_input = [fastq_r2_downsampled_fn]
    else:
        print("Using full FASTQ files (no downsampling).")
        r1_input = [fastq_r1_fn]
        r2_input = [fastq_r2_fn]

    #bc:umi:seq

    print (parse_seqspec_string(args.chemistry))
    protospacer_start, protospacer_len, umi_start, umi_len, barcode_start, barcode_len = parse_seqspec_string(args.chemistry)

    
    # --- 3. Run CRISPR-Correct Mapping ---
    print("Running get_whitelist_reporter_counts_from_fastq...")
    #print (crispr_ambiguous_mapping.__version__, 'version')
    cellbarcode_crisprcorrect_results = crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_fastq(
        whitelist_guide_reporter_df=guide_whitelist_input_df,
        fastq_r1_fns =r1_input ,
        fastq_r2_fns = r2_input,
        # Protospacer parameters (R2)
        protospacer_start_position=protospacer_start,
        protospacer_length=protospacer_len,
        is_protospacer_r1=False,
        is_protospacer_header=False,
        revcomp_protospacer=True,
        protospacer_hamming_threshold_strict=args.hamming_threshold,

        # Guide UMI parameters (R1)
        guide_umi_start_position=umi_start,
        guide_umi_length=umi_len,
        is_guide_umi_r1=True,
        is_guide_umi_header=False,
        revcomp_guide_umi=False,

        # Sample Barcode parameters (R1)
        sample_barcode_start_position=barcode_start,
        sample_barcode_length=barcode_len,
        is_sample_barcode_r1=True,
        is_sample_barcode_header=False,
        revcomp_sample_barcode=False,
        
        cores=args.cores
    )
    print ('starting saving')
    # --- 4. Save Results (Example) ---
    cell_barcode_result = cellbarcode_crisprcorrect_results.all_match_set_whitelist_reporter_counter_series_results.protospacer_match.ambiguous_spread_umi_collapsed_counterseries

    anndata_test = series_to_anndata(cell_barcode_result)
    os.makedirs(f'{args.output_prefix}/counts_unfiltered/', exist_ok=True)
    anndata_test.write_h5ad(f'{args.output_prefix}/counts_unfiltered/adata.h5ad')
    # Clean up downsampled files
    if args.downsample_reads > 0:
        os.remove(fastq_r1_downsampled_fn)
        os.remove(fastq_r2_downsampled_fn)
        print("Cleaned up temporary downsampled files.")

def main():
    print ('its the main')
    parser = argparse.ArgumentParser(
        description="CRISPR guide mapping using crispr_ambiguous_mapping (CRISPR-Correct).",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Input Files ---
    parser.add_argument("-b", "--barcode_inclusion_list_fn", type=str, required=True,
                        default="barcode_inclusion_list_OB1_OB2_merged (1).txt",
                        help="Path to the cell barcode inclusion list file.")
    parser.add_argument("-f", "--fastq", type=str, required=True,
                        default="fastq R1 and R2 separated by space",
                        help="Path to the R1 and R2 FASTQ file (gzip supported). They will be splited")

    parser.add_argument("-g", "--guide_set_fn", type=str, required=True,
                        default="vanilla_guides_18loci_ABE8e_PerturbSeq_guide_set (1)(in).csv",
                        help="Path to the guide set CSV file containing 'spacer' column.")
    parser.add_argument("-o", "--output_prefix", type=str, default="crispr_map_output",
                        help="Prefix for the output result files.")
    
    # --- Processing Parameters ---
    parser.add_argument("-c", "--cores", type=int, default=10,
                        help="Number of CPU cores to use for mapping.")
    parser.add_argument("-d", "--downsample_reads", type=int, default=500000,
                        help="Number of read pairs to downsample for quick testing (set to 0 for full run).")
    
    # --- Mapping Coordinates ---
   
    parser.add_argument("--chemistry", type=str, default="seqspec guide chemistry")

    parser.add_argument("--hamming_threshold", type=int, default=6,
                        help="Maximum Hamming distance for protospacer mapping.")
    
    args = parser.parse_args()
    
    # Replace the default values in the help text with the actual defaults from the original script
    # This block is for display purposes if the user runs with -h, but we can't fully control
    # the help output without modifying the parser constructor, so we'll rely on the defaults above.
    
    run_crispr_mapping(args)

if __name__ == "__main__":
    main()


