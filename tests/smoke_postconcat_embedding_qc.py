"""Container smoke test; no host scientific packages needed."""
import sys
import subprocess
import tempfile
from pathlib import Path
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp
import json

root = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(17)
with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    markers = [line.strip() for line in (root/'assets/cell_cycle/regev_lab_cell_cycle_genes.txt').read_text().splitlines()]
    symbols = markers + [f'CONTROL_{i}' for i in range(30)]
    counts = sp.csr_matrix(rng.poisson(2, size=(90, len(symbols))))
    obs = pd.DataFrame({'batch': ['a']*45+['b']*45, 'pct_counts_mt': [30]*10+[10]*80},
                       index=[f'cell_{i}' for i in range(90)])
    gene = ad.AnnData(counts.copy(), obs=obs.copy(), var=pd.DataFrame({'symbol': symbols}, index=[f'gene_{i}' for i in range(len(symbols))]))
    guide = ad.AnnData(sp.csr_matrix(np.ones((90, 3))), obs=obs.copy())
    guide.layers['guide_assignment'] = sp.csr_matrix(np.tile([1, 0, 0], (90, 1)))
    original = md.MuData({'gene': gene, 'guide': guide})
    original.write_h5mu(d/'input.h5mu')
    subprocess.run([sys.executable, str(root/'bin/postconcat_embedding_qc.py'), str(d/'input.h5mu'),
                    '--output', str(d/'out.h5mu'), '--outdir', str(d/'qc'), '--min-cells-fraction', '0',
                    '--cell-cycle-genes', str(root/'assets/cell_cycle/regev_lab_cell_cycle_genes.txt')], check=True)
    result = md.read_h5mu(d/'out.h5mu')
    assert result.n_obs == 80
    assert (result.mod['gene'].X != counts[10:]).nnz == 0
    assert 'X_pca' not in result.mod['gene'].obsm
    assert 'X_umap' not in result.mod['gene'].obsm
    assert 'log1p' not in result.mod['gene'].uns
    metrics = json.loads((d/'qc/embedding_qc_metrics.json').read_text())
    assert metrics['effective_pcs'] == 50
    assert metrics['effective_neighbors'] == 15
    assert len(list((d/'qc').glob('*.png'))) >= 11
    assert (d/'qc/pca_by_measurement_set.png').exists()
    assert (d/'qc/pca_measurement_sets_colored.png').exists()
    assert (d/'qc/leiden_sweep_umap.png').exists()
    assert metrics['cell_cycle']['status'] == 'completed'
    assert metrics['leiden']['status'] == 'completed'
    print('PASS: MT, PCA facets, Leiden sweep, raw counts and no embeddings in MuData')
