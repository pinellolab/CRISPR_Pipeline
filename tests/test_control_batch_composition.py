"""The control-cell batch composition table the controls evaluation writes.

Motivated by the TAP-seq chr8 screen, where every non-targeting guide had been
delivered in one of fourteen sequencing lanes: 2,033 of 2,049 control-only
cells sat in a single ``obs["batch"]`` level, and the pooled non-targeting
rejection rate was read as a method miscalibration rather than as a batch
confound. These cases check that the table says where the control cells are,
and that the confinement warning fires only when they really are confined.
"""

import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from evaluate_controls import (  # noqa: E402
    COMPOSITION_NOTE,
    COMPOSITION_TABLE,
    control_batch_composition,
    write_control_batch_composition,
)


def _mudata(batches, assignment, targeting, layer=True):
    """A tiny MuData: one gene, three guides, one cell per ``batches`` entry."""
    n_cells = len(batches)
    obs = pd.DataFrame(
        {"batch": pd.Categorical(batches)},
        index=[f"cell{i}" for i in range(n_cells)],
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": [f"guide{i}" for i in range(len(targeting))],
            "intended_target_name": [
                "GENE1" if flag else "non-targeting" for flag in targeting
            ],
            "targeting": list(targeting),
        },
        index=[f"guide{i}" for i in range(len(targeting))],
    )
    matrix = np.asarray(assignment, dtype=np.float32)
    guide = ad.AnnData(
        X=np.zeros_like(matrix) if layer else matrix,
        obs=obs.copy(),
        var=guide_var,
    )
    if layer:
        guide.layers["guide_assignment"] = matrix
    gene = ad.AnnData(
        X=np.ones((n_cells, 1), dtype=np.float32),
        obs=obs.copy(),
        var=pd.DataFrame(index=["GENE1"]),
    )
    return mu.MuData({"gene": gene, "guide": guide})


def test_controls_confined_to_one_batch_are_tabulated_and_warned(tmp_path):
    # lane1 is the small lane that holds every control cell; lane2 is the rest
    # of the screen and carries only targeting guides.
    #   lane1: 2 control-only cells, 1 targeting cell   -> 3 cells
    #   lane2: 0 control-only cells, 7 targeting cells  -> 7 cells
    batches = ["lane1"] * 3 + ["lane2"] * 7
    # guide0/guide1 are non-targeting, guide2 targets GENE1.
    targeting = [False, False, True]
    assignment = [
        [1, 0, 0],  # lane1, control only
        [0, 1, 0],  # lane1, control only
        [0, 0, 1],  # lane1, targeting
    ] + [[0, 0, 1]] * 7  # lane2, targeting

    frame, note, confined = control_batch_composition(
        _mudata(batches, assignment, targeting)
    )

    assert confined is True
    assert list(frame["batch"]) == ["lane1", "lane2"]
    assert list(frame["n_cells"]) == [3, 7]
    assert list(frame["n_control_only_cells"]) == [2, 0]
    # The shares are rounded to six decimals for the report.
    np.testing.assert_allclose(
        frame["control_share_of_batch"], [2 / 3, 0.0], atol=1e-6
    )
    np.testing.assert_allclose(frame["share_of_all_control_cells"], [1.0, 0.0])

    assert "WARNING" in note
    assert "confounded with batch" in note
    assert "'lane1'" in note
    assert "effectively a single batch" in note
    # The note carries both shares the reader needs to see the confound.
    assert "100.0% of all control cells" in note
    assert "30.0% of all cells" in note
    assert "2 level(s)" in note


def test_controls_spread_over_batches_produce_no_warning(tmp_path):
    #   lane1: 2 control-only of 5 cells
    #   lane2: 2 control-only of 5 cells
    batches = ["lane1"] * 5 + ["lane2"] * 5
    targeting = [False, False, True]
    assignment = (
        [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1], [0, 0, 1]]
        + [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1], [0, 0, 1]]
    )

    frame, note, confined = control_batch_composition(
        _mudata(batches, assignment, targeting)
    )

    assert confined is False
    assert sorted(frame["batch"]) == ["lane1", "lane2"]
    assert list(frame["n_cells"]) == [5, 5]
    assert list(frame["n_control_only_cells"]) == [2, 2]
    np.testing.assert_allclose(frame["control_share_of_batch"], [0.4, 0.4])
    np.testing.assert_allclose(frame["share_of_all_control_cells"], [0.5, 0.5])
    assert "WARNING" not in note
    assert "4 cell(s) carry only non-targeting guides" in note


def test_a_cell_with_both_guide_kinds_is_not_a_control_cell():
    batches = ["lane1", "lane1", "lane2"]
    targeting = [False, False, True]
    assignment = [
        [1, 0, 1],  # non-targeting and targeting: not control-only
        [1, 1, 0],  # two non-targeting guides: control-only
        [0, 0, 0],  # no guide at all: not control-only
    ]

    frame, _, confined = control_batch_composition(
        _mudata(batches, assignment, targeting)
    )

    assert confined is False  # lane1 holds all controls but also 2/3 of cells
    assert dict(zip(frame["batch"], frame["n_control_only_cells"])) == {
        "lane1": 1,
        "lane2": 0,
    }


def test_composition_falls_back_to_guide_x_when_there_is_no_layer():
    batches = ["lane1", "lane2"]
    targeting = [False, True]
    frame, _, _ = control_batch_composition(
        _mudata(batches, [[1, 0], [0, 1]], targeting, layer=False)
    )
    assert list(frame["n_control_only_cells"]) == [1, 0]


@pytest.mark.parametrize("fmt", ["csr", "csc"])
def test_composition_reads_a_sparse_assignment(fmt):
    """The pipeline stores the assignment sparse, which is the real input."""
    batches = ["lane1"] * 2 + ["lane2"] * 4
    assignment = [[1, 0, 0], [0, 1, 0]] + [[0, 0, 1]] * 4
    mdata = _mudata(batches, assignment, [False, False, True])
    matrix = getattr(sparse, f"{fmt}_matrix")(
        np.asarray(assignment, dtype=np.float32)
    )
    mdata["guide"].layers["guide_assignment"] = matrix

    frame, _, confined = control_batch_composition(mdata)

    assert list(frame["n_control_only_cells"]) == [2, 0]
    assert confined is True


def test_missing_batch_column_skips_with_a_note(tmp_path):
    mdata = _mudata(["lane1", "lane2"], [[1, 0], [0, 1]], [False, True])
    for modality in mdata.mod.values():
        del modality.obs["batch"]
    mdata.update()
    # ``update`` leaves behind the column it pulled up before the deletion.
    for name in [c for c in mdata.obs.columns if c.endswith("batch")]:
        del mdata.obs[name]

    frame, note, confined = control_batch_composition(mdata)

    assert confined is False
    assert frame.empty
    assert "skipped" in note
    assert "batch" in note
    assert "WARNING" not in note


def test_the_process_publishes_the_whole_evaluation_directory():
    """The new files reach the report only because the output is the directory.

    ``evaluation_controls`` emits ``path "plots"``, so anything the script
    writes into its outdir is captured and the dashboards copy it wholesale.
    Narrowing that to per-file globs would silently drop this table.
    """
    module = (
        REPO_ROOT / "modules" / "local" / "evaluation_controls" / "main.nf"
    ).read_text(encoding="utf-8")
    assert 'path "plots" , emit: evaluation_controls' in module
    assert "*.png" not in module and "*.txt" not in module


def test_write_control_batch_composition_writes_both_files(tmp_path):
    batches = ["lane1"] * 3 + ["lane2"] * 7
    assignment = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] + [[0, 0, 1]] * 7
    mdata = _mudata(batches, assignment, [False, False, True])

    frame, note, confined = write_control_batch_composition(mdata, tmp_path)

    table_path = tmp_path / COMPOSITION_TABLE
    note_path = tmp_path / COMPOSITION_NOTE
    assert table_path.exists() and note_path.exists()

    written = pd.read_csv(table_path, sep="\t")
    assert list(written.columns) == [
        "batch",
        "n_cells",
        "n_control_only_cells",
        "control_share_of_batch",
        "share_of_all_control_cells",
    ]
    pd.testing.assert_frame_equal(
        written.astype({"batch": str}), frame.astype({"batch": str})
    )
    assert note_path.read_text(encoding="utf-8") == note
    assert confined is True
