#!/usr/bin/env python
import argparse

import muon as mu
import numpy as np
import pandas as pd
from gtfparse import read_gtf

from intended_target_key_utils import annotate_intended_target_groups


def _close_mudata_handle(mdata):
    """Close backed MuData file handles when present."""
    file_obj = getattr(mdata, "file", None)
    if file_obj is None:
        return
    close_fn = getattr(file_obj, "close", None)
    if callable(close_fn):
        try:
            close_fn()
        except Exception:
            pass


def generate_pairs_to_test_df(limit, mudata_path, input_gtf, use_backed=True):
    backed_mode = "r" if use_backed else None
    mudata = mu.read_h5mu(mudata_path, backed=backed_mode)
    try:
        guide_data = annotate_intended_target_groups(mudata.mod["guide"].var.copy())
        gene_var_names = pd.Index(mudata.mod["gene"].var_names).astype(str).tolist()
    finally:
        _close_mudata_handle(mudata)

    # Load GTF
    df_gtf_refseq = read_gtf(input_gtf).to_pandas()
    df_gtf_unique_gene_name = df_gtf_refseq.drop_duplicates("gene_name")
    df_gtf_unique_gene_name.set_index("gene_name", inplace=True)

    final_candidates_list = []
    for (
        guide_id,
        target_name,
        target_chr,
        target_start,
        target_end,
        target_key,
        guide_chr,
        guide_start,
    ) in guide_data[
        [
            "guide_id",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
            "intended_target_key",
            "guide_chr",
            "guide_start",
        ]
    ].values:
        if limit is None:
            gene_df = pd.DataFrame(
                {
                    "guide_id": guide_id,
                    "gene_id": gene_var_names,
                    "intended_target_name": target_name,
                    "intended_target_chr": target_chr,
                    "intended_target_start": target_start,
                    "intended_target_end": target_end,
                    "intended_target_key": target_key,
                    "pair_type": "discovery",
                }
            )
            final_candidates_list.append(gene_df)
        else:
            temp_gtf = df_gtf_unique_gene_name.query(f'seqname == "{guide_chr}"').copy()
            temp_gtf["distance_from_guide"] = np.abs(temp_gtf["start"] - guide_start)
            final_candidates = temp_gtf.query(f"distance_from_guide <= {limit}").copy()
            final_candidates["guide_id"] = guide_id
            final_candidates["intended_target_name"] = target_name
            final_candidates["intended_target_chr"] = target_chr
            final_candidates["intended_target_start"] = target_start
            final_candidates["intended_target_end"] = target_end
            final_candidates["intended_target_key"] = target_key
            final_candidates["pair_type"] = "discovery"
            final_candidates.reset_index(inplace=True)
            final_candidates_list.append(
                final_candidates[
                    [
                        "guide_id",
                        "gene_id",
                        "intended_target_name",
                        "intended_target_chr",
                        "intended_target_start",
                        "intended_target_end",
                        "intended_target_key",
                        "pair_type",
                    ]
                ].copy()
            )

    final_candidates_df = pd.concat(final_candidates_list, ignore_index=True)
    final_candidates_df["gene_id"] = final_candidates_df["gene_id"].astype(str).str.split(".").str[0]
    final_candidates_df = final_candidates_df.rename(columns={"gene_id": "gene_name"})
    return final_candidates_df[
        [
            "guide_id",
            "gene_name",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
            "intended_target_key",
            "pair_type",
        ]
    ]


def main(limit, mudata_path, input_gtf, output_path, use_backed=True):
    pairs_df = generate_pairs_to_test_df(
        limit=limit,
        mudata_path=mudata_path,
        input_gtf=input_gtf,
        use_backed=use_backed,
    )
    pairs_df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GTF and guide data.")
    parser.add_argument(
        "--limit",
        type=int,
        default=1000000,
        help="Limit for distance from guide. Use -1 for no limit.",
    )
    parser.add_argument("mudata", type=str, help="Path to the input MuData file")
    parser.add_argument("input_gtf", type=str, help="Path to the input GTF file")
    parser.add_argument(
        "--output",
        type=str,
        default="pairs_to_test.csv",
        help="Output path for generated guide-gene candidate pairs.",
    )
    parser.add_argument(
        "--no_backed",
        action="store_true",
        help="Disable backed MuData reads and load MuData fully in memory.",
    )

    args = parser.parse_args()
    limit = None if args.limit == -1 else args.limit
    main(
        limit=limit,
        mudata_path=args.mudata,
        input_gtf=args.input_gtf,
        output_path=args.output,
        use_backed=not args.no_backed,
    )
