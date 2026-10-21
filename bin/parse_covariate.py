#!/usr/bin/env python

"""Resolve deterministic batch keys and optional samplesheet cov_* columns."""
import json
import re
import sys

import pandas as pd


def add_samplesheet_covariates(batch_table, samplesheet_path):
    sheet = pd.read_csv(samplesheet_path, sep=None, engine="python", dtype=str, keep_default_na=False)
    if "measurement_sets" not in sheet:
        raise ValueError("Samplesheet is missing measurement_sets")
    covariates = [name for name in sheet if name.startswith("cov_")]
    invalid = [name for name in covariates if not re.fullmatch(r"cov_[A-Za-z][A-Za-z0-9_]*", name)]
    if invalid:
        raise ValueError(f"Invalid covariate column name(s): {invalid}")
    if len(set(name.lower() for name in covariates)) != len(covariates):
        raise ValueError("Covariate column names collide ignoring case")
    if not covariates:
        return batch_table
    for name in covariates:
        sheet[name] = sheet[name].str.strip()
        if sheet[name].eq("").any():
            sets = sorted(sheet.loc[sheet[name].eq(""), "measurement_sets"].unique())
            raise ValueError(f"{name} has missing values in measurement set(s): {sets}")
        counts = sheet.groupby("measurement_sets", sort=True)[name].nunique(dropna=False)
        conflicts = sorted(counts.index[counts > 1].tolist())
        if conflicts:
            raise ValueError(f"{name} differs between rows of measurement set(s): {conflicts}")
    metadata = sheet[["measurement_sets", *covariates]].drop_duplicates("measurement_sets")
    if not set(batch_table["batch"].astype(str)).issubset(set(metadata["measurement_sets"])):
        raise ValueError("Prepared batches are absent from the samplesheet covariate table")
    joined = batch_table.merge(metadata, left_on="batch", right_on="measurement_sets", how="left", validate="one_to_one")
    return joined.drop(columns="measurement_sets")


def main():
    batch_table = pd.DataFrame(json.loads(sys.argv[1]))
    if len(sys.argv) > 2:
        batch_table = add_samplesheet_covariates(batch_table, sys.argv[2])
    batch_table.to_csv("parse_covariate.csv", index=False)


if __name__ == "__main__":
    main()
