#!/usr/bin/env python

import argparse
import pandas as pd

def prepare_formula(file_path):
    cov_df = pd.read_csv(file_path)
    final_list = []

    for column in cov_df.columns:
        if len(set(cov_df[column])) > 1:
            final_list.append(column)

    cov_list = ' + '.join(final_list) if final_list else ''

    return cov_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare a formula based on the uniqueness of columns in a CSV file.")
    parser.add_argument("input_file", help="Path to the input parsed covariate CSV file")
    args = parser.parse_args()

    formula = prepare_formula(args.input_file)
    print(formula)
    with open("cov_string.txt", 'w') as f:
        f.write(formula)
