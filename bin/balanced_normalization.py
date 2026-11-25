#!/usr/bin/env python3
"""Visualization-only normalization of balanced raw counts; never filter/write them."""
import argparse
import json
from pathlib import Path
from qc_mudata_io import read_mudata_without_uns
from postconcat_embedding_qc import embedding_plots, qc_observations


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', type=Path, default=Path('embedding_qc'))
    p.add_argument('--batch-key', default='batch')
    p.add_argument('--seed', type=int, default=20261026)
    p.add_argument('--cell-cycle-genes', type=Path, required=True)
    args = p.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    args.stage='balanced_raw_counts'; args.n_pcs=50; args.n_neighbors=15; args.n_top_genes=6000
    args.hvg_batch_key=''; args.max_dense_gb=24; args.cell_cycle='auto'; args.reference='human'
    args.tapseq_mode=False; args.min_cell_cycle_genes=10; args.enable_leiden=True
    args.leiden_resolutions=[0.2,0.5,1.0]; args.leiden_diagnostic_resolution=0.5; args.leiden_n_iterations=2
    data=read_mudata_without_uns(args.input)
    report=embedding_plots(data['gene'],qc_observations(data,args.batch_key),args)
    report.update(stage=args.stage, input_cells=data.n_obs, retained_cells=data.n_obs,
                  normalized_matrix_saved=False, cell_filtering=False, gene_filtering=False,
                  scope='Temporary visualization copy; raw balanced H5MU is unchanged')
    (args.outdir/'embedding_qc_metrics.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__': main()
