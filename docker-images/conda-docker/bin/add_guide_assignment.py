#!/usr/bin/env python
import muon as mu
import argparse
from scipy.io import mmread
from scipy.sparse import csr_matrix
import numpy as np

def add_guide_assignment(mudata_path, guide_assignment_mtx):
    # Load MuData object
    mudata = mu.read_h5mu(mudata_path)
    
    sparse_matrix = mmread(guide_assignment_mtx).T
    sparse_matrix_csr = csr_matrix(sparse_matrix)  # Convert to CSR format
    
    # Add to mudata
    mudata.mod['guide'].layers['guide_assignment'] = sparse_matrix_csr
    
    # Save MuData
    mudata.write("sceptre_assignment_mudata.h5mu")
    print(f"Successfully added guide_assignment to mudata and saved as sceptre_assignment_mudata.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process MuData and guide assignment.')
    parser.add_argument('--mudata', required=True, help='Path to the input h5mu file.')
    parser.add_argument('--guide_assignment', required=True, 
                        help='Path to the guide assignment in Matrix Market format (.mtx).')
    
    args = parser.parse_args()
    add_guide_assignment(args.mudata, args.guide_assignment)