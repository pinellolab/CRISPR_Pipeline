import importlib.util
from pathlib import Path

import numpy as np
import scipy.sparse as sp


SCRIPT = Path(__file__).parents[1] / "bin" / "remove_clonal_cells.py"
SPEC = importlib.util.spec_from_file_location("remove_clonal_cells", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_published_grouping_and_filter_actions():
    rows = [
        list(range(0, 10)),
        list(range(0, 10)),
        list(range(10, 20)),
        list(range(0, 20)),
    ]
    row_index = np.repeat(np.arange(4), [len(row) for row in rows])
    col_index = np.concatenate(rows)
    assignment = sp.csr_matrix(
        (np.ones(len(col_index)), (row_index, col_index)), shape=(4, 1000)
    )
    grouped = MODULE.group_clones(assignment, alpha=0.05)

    assert grouped["members"] == [[0, 1], [2]]
    assert grouped["ambiguous"].tolist() == [False, False, False, True]

    keep, clonal, sizes, representative = MODULE.decide_cells(
        grouped,
        action="drop_clonal",
        min_clone_size=2,
        gene_totals=np.array([100, 200, 150, 300]),
    )
    assert keep.tolist() == [False, False, True, False]
    assert clonal.tolist() == [True, True, False, False]
    assert sizes.tolist() == [2, 2, 1, 1]
    assert representative.tolist() == [False, True, True, False]

    keep_one, *_ = MODULE.decide_cells(
        grouped,
        action="keep_representative",
        min_clone_size=2,
        gene_totals=np.array([100, 200, 150, 300]),
    )
    assert keep_one.tolist() == [False, True, True, False]

