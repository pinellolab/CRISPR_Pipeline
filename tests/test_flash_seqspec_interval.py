"""Shifted spacers must not steal a feature extracted at a known position."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import flash_base_editing_mapping as flash


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_seqspec_interval_recovers_shifted_spacer_identity(backend, monkeypatch):
    if backend == 'numpy':
        monkeypatch.setattr(flash, 'torch', None)
    elif flash.torch is None:
        pytest.skip('torch is not installed')
    guides = ['AAACCC', 'AACCCG']
    matcher = flash.DenseGuideMatcher(guides, 6, 0, 'cpu', 32, unique_only=True)
    read = 'AAACCCG'
    # Legacy sliding matching sees both exact sequences; interval 1:7 identifies guide 2.
    extracted = flash.slice_from_read_pair('', read, (1, 1, 7))
    _, start, gid = matcher.best_hits([extracted])
    assert gid.tolist() == [1]
    assert start.tolist() == [0]


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_ambiguous_mismatch_is_unassigned_and_exact_hit_wins(backend, monkeypatch):
    if backend == 'numpy':
        monkeypatch.setattr(flash, 'torch', None)
    elif flash.torch is None:
        pytest.skip('torch is not installed')
    matcher = flash.DenseGuideMatcher(['AAAA', 'AACC'], 4, 1, 'cpu', 32, unique_only=True)
    distances, _, gids = matcher.best_hits(['AAAC', 'AAAA', 'AA', 'TTTT'])
    assert gids.tolist() == [-1, 0, -1, -1]
    assert distances[1] == 0
    assert matcher.ambiguous_reads == 1


def test_legacy_full_read_tie_policy_is_preserved(monkeypatch):
    monkeypatch.setattr(flash, 'torch', None)
    matcher = flash.DenseGuideMatcher(['AAAA', 'AACC'], 4, 1, 'cpu', 32)
    _, _, gids = matcher.best_hits(['AAAC'])
    assert gids.tolist() == [0]
