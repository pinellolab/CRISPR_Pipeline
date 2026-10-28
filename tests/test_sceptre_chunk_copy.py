"""A byte-identical calibration handoff and failure-safe publication."""
import sys
from pathlib import Path

import anndata as ad
import h5py
import mudata as md
import numpy as np
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import chunk_mudata_sceptre as chunker


def test_single_chunk_preserves_full_hdf5(tmp_path):
    cells = [f'cell{i}' for i in range(6)]
    gene = ad.AnnData(sparse.csr_matrix(np.arange(24).reshape(6, 4)))
    guide = ad.AnnData(sparse.csr_matrix(np.ones((6, 2))))
    gene.obs_names = guide.obs_names = cells
    guide.layers['guide_assignment'] = guide.X.copy()
    data = md.MuData({'gene': gene, 'guide': guide})
    data.uns['balanced_matching'] = {'seed': 42, 'rna_umis': int(gene.X.sum())}
    source = tmp_path/'source.h5mu'; data.write_h5mu(source)
    files, _, summary = chunker.chunk_mudata_sceptre(str(source), str(tmp_path/'chunks'), chunk_mode='off')
    assert not summary['chunked']
    assert source.read_bytes() == Path(files[0]).read_bytes()
    with h5py.File(files[0]) as f:
        f.visititems(lambda name, node: None)
    restored = md.read_h5mu(files[0])
    assert restored.uns['balanced_matching']['rna_umis'] == gene.X.sum()
    assert (restored['guide'].layers['guide_assignment'] != guide.X).nnz == 0


def test_checksum_failure_does_not_replace_previous_output(tmp_path, monkeypatch):
    source = tmp_path/'input'; source.write_bytes(b'complete input')
    target = tmp_path/'output'; target.write_bytes(b'previous output')
    monkeypatch.setattr(chunker, '_file_sha256', lambda path: 'invalid checksum')
    with pytest.raises(OSError, match='checksum mismatch'):
        chunker._verified_stream_copy(source, target)
    assert target.read_bytes() == b'previous output'
    assert not list(tmp_path.glob('.chunk-copy-*'))
