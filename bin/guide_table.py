#!/usr/bin/env python

import pandas as pd
import argparse


def reverse_complement(seq):
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

def process_table(guide_table, rev_comp):
    if rev_comp.lower() == 'true':
        rev_comp = True
    elif rev_comp.lower() == 'false':
        rev_comp = False  
    else:
        raise ValueError("rev_comp must be 'true' or 'false'")  



    guide_metadata = pd.read_csv(guide_table, sep="\t")
    if rev_comp:
        print ('using the reverse complement option')
        # Function to get reverse complement
        guide_metadata['spacer'] = guide_metadata['spacer'].apply(reverse_complement)
    
    guide_metadata[['spacer', 'guide_id']].to_csv('guide_features.txt',
                                                           sep='\t', header=None, index=None)

def main():
    parser = argparse.ArgumentParser(description="Process guide metadata.")
    parser.add_argument('--guide_table', type=str, required=True, help="Path to the input Guide file.")
    #boolean flag for reverse complement option
    parser.add_argument('--rev_comp', type=str, required=True,  help="Flag to indicate if reverse complement is needed.")
    args = parser.parse_args()
    process_table(args.guide_table, args.rev_comp) 

if __name__ == "__main__":
    main()
