import pathlib
import sys

import pandas as pd


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from prepare_inference import _pairs_to_h5mu_dict


def test_pairs_to_h5mu_dict_normalizes_mixed_chromosomes():
    pairs = pd.DataFrame(
        {
            "guide_id": ["g1", "g2", "g3"],
            "gene_name": ["A", "B", "C"],
            "intended_target_chr": ["chr8", 8, None],
        }
    )

    observed = _pairs_to_h5mu_dict(pairs)

    assert observed["intended_target_chr"] == ["chr8", "8", ""]
