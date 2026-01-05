#!/usr/bin/env python

import pandas as pd
import argparse


def reverse_complement(seq):
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

def process_table(guide_table, rev_comp, spacer):
    print (spacer, 'spacer tag function ')
    if spacer == '':
        print ('no spacer tag added to the guide sequence')
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

    #adding a tag in the guide sequence
    if spacer:
        print ('spacer tag added to the guide sequence', 'spacer tag:', spacer  )
        guide_metadata['spacer'] = guide_metadata['spacer'].apply(lambda x : spacer + x )

    guide_metadata[['spacer', 'guide_id']].to_csv('guide_features.txt',
                                                           sep='\t', header=None, index=None)
    
    read_guide_features = pd.read_csv('guide_features.txt', sep="\t", header=None)
    print (read_guide_features.head())

def main():
    parser = argparse.ArgumentParser(description="Process guide metadata.")
    parser.add_argument('--guide_table', type=str, required=True, help="Path to the input Guide file.")
    #boolean flag for reverse complement option
    parser.add_argument('--rev_comp', type=str, required=True,  help="Flag to indicate if reverse complement is needed.")
    parser.add_argument('--spacer', type=str, default='', help="Tag to add to the guide sequence. Ex: 'TAGCTCTTAAAC'")
    args = parser.parse_args()
    process_table(args.guide_table, args.rev_comp, args.spacer) 

if __name__ == "__main__":
    main()
