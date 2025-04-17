#!/usr/bin/env python

import argparse
import muon as mu
import anndata as ad
import pandas as pd
import numpy as np
from muon import MuData

import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser(description="Scrublet doublet detection for MuData object")
    parser.add_argument('--input', type=str, required=True, help='Path to the input MuData file')

    args = parser.parse_args()
    mdata = mu.read(args.input)

    # Run scrublet
    scrub = scr.Scrublet(mdata.mod['gene'].X, random_state=42)
    mdata.mod['gene'].obs['doublet_scores'], mdata.mod['gene'].obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.plot_histogram()

    # Save plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")
    plt.savefig(f"figures/doublets_batch.png")
    
    print("Number of predicted doublets:", sum(mdata.mod['gene'].obs['predicted_doublets']))

    mdata.mod['gene'].obs['doublet_info'] = mdata.mod['gene'].obs["predicted_doublets"].astype(str)

    # Remove doublets
    mdata.mod['gene'] = mdata.mod['gene'][mdata.mod['gene'].obs['predicted_doublets'] == False]
    print("Doublets removed.")

    # Update mudata
    # Find the intersection of barcodes between scRNA and guide data
    intersecting_barcodes = list(set(mdata.mod['gene'].obs_names)
                                 .intersection(mdata.mod['guide'].obs_names))

    mdata = MuData({
        'gene': mdata.mod['gene'][intersecting_barcodes, :].copy(),
        'guide': mdata.mod['guide'][intersecting_barcodes, :].copy()
    })

    obs_names = set(mdata.mod['guide'].obs.columns.tolist()) & set(mdata.mod['gene'].obs.columns.tolist())
    mdata.obs = mdata.mod['guide'].obs.loc[:, list(obs_names)]

    # Save
    mdata.write('mdata_doublets.h5mu')

if __name__ == '__main__':
    main()
