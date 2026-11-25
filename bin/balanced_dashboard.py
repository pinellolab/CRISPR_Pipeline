#!/usr/bin/env python3
"""Render a fresh balanced-data scientific report and canonical QC metrics."""
import argparse
import html
import json
from pathlib import Path

import numpy as np
from qc_mudata_io import read_mudata_without_uns
from qc_metrics_json import build_pipeline_qc_metrics_payload
from dashboard_presentation import enhance_dashboard_html
from balance_mudata import signature


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--params', type=Path, required=True)
    args = p.parse_args()
    report = json.loads(args.report.read_text())
    params = json.loads(args.params.read_text())
    data = read_mudata_without_uns(args.input)
    if (data.n_obs != report['selected_cells']
            or int(data['gene'].X.sum(dtype=np.int64)) != report['analysis_gene_rna_umis']
            or signature(data['guide'].X) != report['guide_counts_sha256']
            or signature(data['guide'].layers['guide_assignment']) != report['guide_assignment_sha256']):
        raise ValueError('Final inference object differs from the balanced counts/cells/guide calls')
    n = data.n_obs
    payload = build_pipeline_qc_metrics_payload(
        params=params, qc_params={'params_found': True, 'barcode_filter': 'inherited parent cohort'},
        sample_counts=[], concat_count=n, filter_count=n, filter_info={}, filtered_count=n,
        guide_intersection_count=n, final_count=n, mudata=data,
        gene_ann=data['gene'], gene_filtered_ann=data['gene'], guide_ann=data['guide'],
        additional_qc_dir='additional_qc', use_default=True)
    payload['balanced_matching'] = report
    payload['sources']['dashboard_builder'] = 'balanced_dashboard.py'
    payload['sources']['upstream_qc_scope'] = report['source_qc']
    Path('pipeline_qc_metrics.json').write_text(json.dumps(payload, indent=2) + '\n')
    esc = html.escape
    rows = ''.join(f'<tr><th>{esc(str(key))}</th><td>{esc(str(value))}</td></tr>' for key, value in report.items())
    parts = [f'<html><head><meta charset="utf-8"><title>{esc(report["dataset"])} balanced benchmark</title>'
             '<style>body{font-family:system-ui;max-width:1400px;margin:30px auto;padding:20px}img{max-width:100%}'
             'td,th{padding:8px;text-align:left;overflow-wrap:anywhere}figure{border-top:1px solid #ddd;padding:20px}</style>'
             '</head><body>', f'<h1>{esc(report["dataset"])}: balanced RNA benchmark</h1>',
             '<p>Fresh inference from hypergeometrically thinned raw RNA; guide counts and assignments are fixed. '
             'Parent cohort eligibility is inherited. Mapping, doublet calls and guide assignment were not repeated. '
             'Full-transcriptome matching precedes the frozen parent gene subset used for inference. '
             'Median-depth normalization, log1p, HVGs, scaling, PCA and UMAP affect a temporary visualization copy only. '
             'Equal cell counts and global UMI totals do not remove biological or assay differences.</p>',
             f'<table>{rows}</table>']
    for folder in ('additional_qc', 'evaluation_output', 'control_output', 'benchmark_output'):
        parts.append(f'<h2>{esc(folder.replace("_", " "))}</h2>')
        for path in sorted(Path(folder).rglob('*')):
            if not path.is_file() or path.name.startswith('.'): continue
            name = path.as_posix(); suffix = path.suffix.lower()
            if suffix in ('.png', '.svg', '.jpg', '.jpeg'):
                parts.append(f'<figure><figcaption>{esc(name)}</figcaption><img src="{esc(name, quote=True)}" alt="{esc(name, quote=True)}"></figure>')
            elif suffix in ('.json', '.tsv', '.csv', '.html', '.pdf'):
                parts.append(f'<p><a href="{esc(name, quote=True)}">{esc(name)}</a></p>')
    parts.append('</body></html>')
    Path('dashboard.html').write_text(enhance_dashboard_html(''.join(parts), run_status='completed'))


if __name__ == '__main__': main()
