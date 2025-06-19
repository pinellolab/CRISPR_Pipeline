#!/usr/bin/env python
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', default='updated_samplesheet.tsv')
args = parser.parse_args()

samplesheet = pd.read_csv(args.input, sep='\t')

updated_samplesheet = samplesheet[['R1_path', 'R2_path', 'file_modality', 'measurement_sets', 'sequencing_run', 'lane', 'flowcell_id','seqspec', 'barcode_onlist']].copy()

updated_samplesheet['file_modality'] = updated_samplesheet['file_modality'].replace({
    'scRNA sequencing': 'scRNA',
    'gRNA sequencing': 'gRNA',
    'cell hashing barcode sequencing': 'hash'
})

updated_samplesheet.to_csv(args.output, sep=',', index=False)
