#!/usr/bin/env python
import mudata as mu
import argparse
from scipy.io import mmread
import numpy as np

from count_matrix_utils import describe_matrix, to_sparse_counts

def add_guide_assignment(mudata_path, guide_assignment_mtx):
    # Load MuData object
    mudata = mu.read_h5mu(mudata_path)

    sparse_matrix = mmread(guide_assignment_mtx).T
    # Matrix Market "real" headers make mmread return float64 even for the
    # 0/1 assignment calls this file holds; store CSR at the narrowest dtype
    # that fits the data instead (uint16 in practice).
    sparse_matrix_csr = to_sparse_counts(sparse_matrix)
    print(describe_matrix(sparse_matrix_csr, "guide_assignment"))

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