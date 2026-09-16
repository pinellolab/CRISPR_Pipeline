import pathlib
import subprocess
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

import evaluate_controls as ec
from evaluate_controls import (
    perform_binary_evaluation,
    run_evaluation_controls,
    run_evaluation_controls_from_path,
)


def test_control_evaluation_skips_cleanly_without_negative_controls(tmp_path):
    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(
        X=np.ones((1, 1)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["GENE1"]),
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": ["guide1"],
            "intended_target_name": ["GENE1"],
            "targeting": [True],
        },
        index=["guide1"],
    )
    guide = ad.AnnData(X=np.ones((1, 1)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["guide1"],
            "gene_id": ["GENE1"],
            "perturbo_log2_fc": [-1.0],
            "perturbo_p_value": [0.01],
        }
    )

    run_evaluation_controls(mdata, outdir=tmp_path)

    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "no non-targeting guides" in marker.read_text()


def test_binary_evaluation_skips_empty_input(tmp_path):
    assert not perform_binary_evaluation([], [], tmp_path, plot=False)
    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "valid_rows=0" in marker.read_text()


def test_control_evaluation_matches_symbol_target_to_ensembl_result(tmp_path, monkeypatch):
    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(X=np.ones((1, 1)), obs=obs.copy(), var=pd.DataFrame(index=["ENSG1"]))
    guide_var = pd.DataFrame(
        {
            "guide_id": ["targeting-guide", "control-guide"],
            "intended_target_name": ["GENE1", "non-targeting"],
            "targeting": [True, False],
        },
        index=["targeting-guide", "control-guide"],
    )
    guide = ad.AnnData(X=np.ones((1, 2)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["targeting-guide", "control-guide"],
            "gene_id": ["ENSG1", "ENSG1"],
            "gene_name": ["GENE1", "GENE1"],
            "perturbo_log2_fc": [-1.0, 0.0],
            "perturbo_p_value": [0.01, 0.9],
        }
    )
    monkeypatch.setattr("evaluate_controls.plot_volcano", lambda *args, **kwargs: None)
    monkeypatch.setattr("evaluate_controls.perform_binary_evaluation", lambda *args, **kwargs: True)
    monkeypatch.setattr("evaluate_controls.savefig", lambda *args, **kwargs: None)

    run_evaluation_controls(mdata, outdir=tmp_path)

    assert not (tmp_path / "controls_evaluation_skipped.txt").exists()


def test_untested_pairs_do_not_unbalance_the_matched_classes(tmp_path, monkeypatch):
    import matplotlib

    matplotlib.use("Agg", force=True)

    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(
        X=np.ones((1, 2)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["ENSG1", "ENSG2"]),
    )
    guide_ids = ["t1", "t2", "c1", "c2", "c3"]
    guide_var = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "intended_target_name": [
                "GENE1",
                "GENE2",
                "non-targeting|1",
                "non-targeting|1",
                "non-targeting|2",
            ],
            "targeting": [True, True, False, False, False],
        },
        index=guide_ids,
    )
    guide = ad.AnnData(X=np.ones((1, 5)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    # t2/ENSG2 is a direct target the conditional randomization test did not
    # test, so it carries no p-value; c3/ENSG1 is an untested control.
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["t1", "t2", "c1", "c2", "c3"],
            "gene_id": ["ENSG1", "ENSG2", "ENSG1", "ENSG2", "ENSG1"],
            "gene_name": ["GENE1", "GENE2", "GENE1", "GENE2", "GENE1"],
            "perturbo_log2_fc": [-1.0, -1.0, 0.0, 0.0, 0.0],
            "perturbo_p_value": [0.01, np.nan, 0.4, 0.5, np.nan],
        }
    )
    monkeypatch.setattr("evaluate_controls.plot_volcano", lambda *args, **kwargs: None)
    monkeypatch.setattr("evaluate_controls.savefig", lambda *args, **kwargs: None)

    run_evaluation_controls(mdata, outdir=tmp_path)

    summary = (tmp_path / "controls_evaluation_summary.txt").read_text()
    assert "direct_target_rows=1" in summary
    assert "non_targeting_rows=1" in summary
    # The classes are matched after the untested pairs are removed, so nothing
    # is left for the curve builder to drop.
    assert "unscorable_rows_dropped=0" in summary
    assert "positive_rate=0.500000" in summary


def test_control_evaluation_skips_cleanly_without_global_results(tmp_path):
    class LocalOnlyResult:
        uns = {"local_analysis_per_guide_results": pd.DataFrame()}

    run_evaluation_controls(LocalOnlyResult(), outdir=tmp_path)

    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "global PerTurbo guide results are not present" in marker.read_text()


# --------------------------------------------------------------------------
# The evaluation reads the result table out of the file rather than loading it,
# and reads most of its columns only for the rows that reach a plot. What that
# must not change is which rows those are, so the cases below are pinned to the
# files the previous whole-frame implementation emitted for the same input
# (bin/evaluate_controls.py at 7992398, run side by side over these fixtures).
# --------------------------------------------------------------------------

GUIDES = [
    # guide_id, intended_target_name, targeting
    ("t1", "GENE1", True),
    ("t2", "ensg2.7", True),      # lowercase, with an Ensembl version suffix
    ("t3", " ENSG3 ", True),      # whitespace around the identifier
    ("t4", "", True),             # blank target: carries no identifier
    ("c1", "non-targeting", False),
    ("c2", "non-targeting", False),
    ("c3", None, False),          # no intended target at all
]

# The same library aimed at genomic elements instead of genes, as an enhancer
# screen is: nothing a guide intends is a tested gene, so no row is a direct
# target and the evaluation has nothing to score.
GUIDES_ELEMENTS = [
    ("t1", "chr8:100-200", True),
    ("t2", "chr8:300-400", True),
    ("t3", "chr8:500-600", True),
    ("t4", "chr8:700-800", True),
    ("c1", "non-targeting", False),
    ("c2", "non-targeting", False),
    ("c3", None, False),
]

# guide_id, gene_id, gene_name, p_value
ROWS = [
    ("t1", "ENSG1", "GENE1", 0.001),
    ("t1", "ENSG2", "GENE2", 0.5),
    ("t2", "ENSG2", "GENE2", 0.002),
    ("t2", "ENSG1", "GENE1", 0.6),
    ("t3", "ENSG3", "GENE3", np.nan),    # a direct target the CRT never tested
    ("t3", "ENSG4", "GENE4", 0.7),
    ("t4", "ENSG4", "GENE4", 0.01),
    ("orphan", "ENSG1", "GENE1", 0.02),  # guide absent from the guide metadata
    ("c1", "ENSG1", "GENE1", 0.4),       # control tested against a direct target
    ("c1", "ENSG4", "GENE4", 0.45),      # control tested against something else
    ("c2", "ENSG2", "GENE2", 0.55),
    ("c2", "ENSG4", "GENE4", np.nan),
    ("c3", "ENSG1", "GENE1", np.nan),    # matched control with no p-value
    ("c3", "ENSG2", "GENE2", 0.65),
]

# No control was ever tested against a direct target, so the evaluation falls
# back to a random sample of the controls.
ROWS_NO_MATCHED_CONTROLS = [row for row in ROWS if not row[0].startswith("c")] + [
    ("c1", "ENSG4", "GENE4", 0.45),
    ("c2", "ENSG4", "GENE4", 0.35),
    ("c3", "ENSG4", "GENE4", 0.25),
]

# More direct targets than matched controls: the controls are sampled with
# replacement.
ROWS_REPLACEMENT = [
    ("t1", "ENSG1", "GENE1", 0.001),
    ("t2", "ENSG2", "GENE2", 0.002),
    ("t3", "ENSG3", "GENE3", 0.003),
    ("c1", "ENSG1", "GENE1", 0.4),
    ("c2", "ENSG2", "GENE2", 0.55),
]

# Classes the p-value does not separate, and more matched controls than direct
# targets, so the areas depend on which controls the sample drew -- these are
# the numbers that move if the control rows are counted or ordered differently.
ROWS_INTERLEAVED = [
    ("t1", "ENSG1", "GENE1", 0.001),
    ("t2", "ENSG2", "GENE2", 0.30),
    ("t3", "ENSG3", "GENE3", 0.60),
    ("c1", "ENSG1", "GENE1", 0.10),
    ("c1", "ENSG2", "GENE2", 0.05),
    ("c2", "ENSG2", "GENE2", 0.50),
    ("c3", "ENSG1", "GENE1", 0.20),
    ("c3", "ENSG2", "GENE2", 0.70),
]

# Two direct-target rows that agree on every column the evaluation reads and
# differ only in one it does not. drop_duplicates ran over the whole frame, so
# both survive, and the scored counts say so.
ROWS_DUPLICATES = [
    ("t1", "ENSG1", "GENE1", 0.001),
    ("t1", "ENSG1", "GENE1", 0.001),
    ("t2", "ENSG2", "GENE2", 0.002),
    ("c1", "ENSG1", "GENE1", 0.4),
    ("c2", "ENSG2", "GENE2", 0.55),
    ("c3", "ENSG1", "GENE1", 0.6),
]

ROWS_NO_P_VALUES = [(g, gi, gn, np.nan) for g, gi, gn, _ in ROWS]

PLOTS = [
    "global_analysis_perturbo_barplot_direct_vs_control.png",
    "global_analysis_perturbo_precision_recall_roc.png",
    "global_analysis_perturbo_volcano_plot.png",
]

SUMMARY = "controls_evaluation_summary.txt"
SKIPPED = "controls_evaluation_skipped.txt"

# Written for every run, whatever the result table holds: where the control
# cells sit across obs["batch"] decides how the calibration numbers read.
COMPOSITION = [
    "control_batch_composition.tsv",
    "control_batch_composition_note.txt",
]


def _summary(direct, non_targeting, auprc, auroc, dropped=0):
    return (
        f"direct_target_rows={direct}\n"
        f"non_targeting_rows={non_targeting}\n"
        f"unscorable_rows_dropped={dropped}\n"
        f"positive_rate={direct / (direct + non_targeting):.6f}\n"
        f"auprc={auprc}\n"
        f"auroc={auroc}\n"
    )


def _skip(group, rows, controls):
    return (
        f"Controls evaluation skipped: no {group} with a finite "
        "perturbo_p_value are present in the inference results, so AUROC/AUPRC "
        "and matched-control plots cannot be calculated.\n"
        f"targeting_direct_target_rows={rows}\n"
        f"non_targeting_control_rows={controls}\n"
    )


# case -> (mudata builder kwargs, expected file, expected content, expected plots)
EXPECTED = {
    "mixed": (
        dict(rows=ROWS),
        SUMMARY,
        _summary(2, 2, "1.000000", "1.000000"),
        PLOTS,
    ),
    # the guide table arrives with "TRUE"/"FALSE" instead of booleans
    "targeting_as_text": (
        dict(rows=ROWS, targeting_as_text=True),
        SUMMARY,
        _summary(2, 2, "1.000000", "1.000000"),
        PLOTS,
    ),
    # identifiers stored as plain strings rather than as categories
    "plain_strings": (
        dict(rows=ROWS, categorical=False),
        SUMMARY,
        _summary(2, 2, "1.000000", "1.000000"),
        PLOTS,
    ),
    "no_matched_controls": (
        dict(rows=ROWS_NO_MATCHED_CONTROLS),
        SUMMARY,
        _summary(2, 2, "1.000000", "1.000000"),
        PLOTS,
    ),
    "replacement": (
        dict(rows=ROWS_REPLACEMENT),
        SUMMARY,
        _summary(3, 3, "1.000000", "1.000000"),
        PLOTS,
    ),
    "interleaved": (
        dict(rows=ROWS_INTERLEAVED),
        SUMMARY,
        _summary(3, 3, "0.711111", "0.666667"),
        PLOTS,
    ),
    "duplicate_rows": (
        dict(rows=ROWS_DUPLICATES, unique_unread=True),
        SUMMARY,
        _summary(3, 3, "1.000000", "1.000000"),
        PLOTS,
    ),
    "no_gene_name": (
        dict(rows=ROWS, drop=("gene_name",)),
        SUMMARY,
        _summary(1, 1, "1.000000", "1.000000"),
        PLOTS,
    ),
    "no_target_column": (
        dict(rows=ROWS, drop=("intended_target_name",)),
        SUMMARY,
        _summary(2, 2, "1.000000", "1.000000"),
        PLOTS,
    ),
    "all_p_missing": (
        dict(rows=ROWS_NO_P_VALUES),
        SKIPPED,
        _skip("direct-target rows", 0, 0),
        [],
    ),
    "no_direct_targets": (
        dict(rows=ROWS, elements=True),
        SKIPPED,
        _skip("direct-target rows", 0, 4),
        [],
    ),
    "no_controls": (
        dict(rows=[r for r in ROWS if not r[0].startswith("c")]),
        SKIPPED,
        _skip("non-targeting guides", 2, 0),
        [],
    ),
    "absent_key": (
        dict(rows=None),
        SKIPPED,
        "Controls evaluation skipped: global PerTurbo guide results are not "
        "present for this inference mode.\n"
        "global_analysis_per_guide_results=absent\n",
        [],
    ),
}


def _results_frame(rows, categorical=True, drop=(), unique_unread=False):
    rng = np.random.default_rng(7)
    n = len(rows)
    frame = pd.DataFrame(
        {
            "gene_id": [r[1] for r in rows],
            "guide_id": [r[0] for r in rows],
            "perturbo_log2_fc": np.asarray(
                [-1.5, 0.2, -2.0, 0.1, -1.1, 0.3, -0.9, 0.4,
                 0.05, -0.1, 0.2, -0.3, 0.15, -0.2][:n],
                dtype=np.float32,
            ),
            "perturbo_p_value": [r[3] for r in rows],
            "perturbo_q_value": rng.random(n),
            # a stale target name the evaluation is expected to overwrite
            "intended_target_name": ["stale"] * n,
            "targeting": [not r[0].startswith("c") for r in rows],
            "gene_name": [r[2] for r in rows],
            "nPerturbedCells": pd.array(
                [10 + i if i % 3 else None for i in range(n)], dtype="Int64"
            ),
            # a column the evaluation never reads
            "unused_metric": np.arange(n, dtype=float) if unique_unread else np.zeros(n),
        }
    )
    if categorical:
        for column in ("gene_id", "guide_id", "gene_name", "intended_target_name"):
            frame[column] = frame[column].astype("category")
    return frame.drop(columns=list(drop))


def _write_fixture(
    tmp_path,
    rows,
    categorical=True,
    drop=(),
    unique_unread=False,
    targeting_as_text=False,
    elements=False,
):
    """An inference MuData shaped like the screens, four genes wide."""
    guides = GUIDES_ELEMENTS if elements else GUIDES
    guide_ids = [g[0] for g in guides]
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(4)])
    gene = ad.AnnData(
        X=sparse.csr_matrix(np.arange(16, dtype=np.float32).reshape(4, 4)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["ENSG1", "ENSG2", "ENSG3", "ENSG4"]),
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "intended_target_name": pd.Categorical([g[1] for g in guides]),
            "targeting": (
                [str(g[2]).upper() for g in guides]
                if targeting_as_text
                else [g[2] for g in guides]
            ),
        },
        index=guide_ids,
    )
    guide = ad.AnnData(
        X=sparse.csr_matrix(np.eye(4, len(guide_ids), dtype=np.float32)),
        obs=obs.copy(),
        var=guide_var,
    )
    mdata = mu.MuData({"gene": gene, "guide": guide})
    if rows is not None:
        mdata.uns["global_analysis_per_guide_results"] = _results_frame(
            rows, categorical=categorical, drop=drop, unique_unread=unique_unread
        )
    # a second result table, to be left unread
    mdata.uns["local_analysis_per_guide_results"] = _results_frame(ROWS[:4])
    path = tmp_path / "inference_mudata.h5mu"
    mdata.write(path)
    return path


@pytest.fixture(autouse=True)
def _headless():
    import matplotlib

    matplotlib.use("Agg", force=True)


@pytest.mark.parametrize("case", list(EXPECTED))
def test_reading_less_emits_the_same_metrics(tmp_path, case):
    kwargs, name, expected, plots = EXPECTED[case]
    path = _write_fixture(tmp_path, **kwargs)
    outdir = tmp_path / "plots"

    run_evaluation_controls_from_path(str(path), str(outdir))

    assert (outdir / name).read_text(encoding="utf-8") == expected
    assert sorted(p.name for p in outdir.iterdir()) == sorted(
        [name, *plots, *COMPOSITION]
    )


def test_the_evaluation_never_loads_the_result_tables(tmp_path, monkeypatch):
    """read_h5mu would read all of uns; on the real screens that is 6.35 GB."""
    path = _write_fixture(tmp_path, rows=ROWS)

    def refuse(*args, **kwargs):
        raise AssertionError("read_h5mu loads every result table in uns")

    monkeypatch.setattr(mu, "read_h5mu", refuse)

    run_evaluation_controls_from_path(str(path), str(tmp_path / "plots"))

    assert (tmp_path / "plots" / SUMMARY).exists()


def test_a_large_row_selection_reads_the_same_values(tmp_path, monkeypatch):
    """Past a point, columns are read in order and masked rather than by point.

    Real screens cross that threshold -- a screen whose controls were all
    tested against the intended targets selects millions of rows -- so the two
    readers have to return the same thing.
    """
    monkeypatch.setattr(ec, "_POINT_SELECT_LIMIT", 1)
    monkeypatch.setattr(ec, "_SLAB_ROWS", 3)
    path = _write_fixture(tmp_path, rows=ROWS_INTERLEAVED)
    outdir = tmp_path / "plots"

    run_evaluation_controls_from_path(str(path), str(outdir))

    _, name, expected, _ = EXPECTED["interleaved"]
    assert (outdir / name).read_text(encoding="utf-8") == expected


def test_duplicate_rows_differing_only_in_an_unread_column_are_kept(tmp_path):
    """drop_duplicates ran over every column, so a difference anywhere counts.

    Both t1/ENSG1 rows agree on the guide, the gene, the fold change and the
    p-value and differ only in a column nothing reads, and the evaluation
    scores three direct targets rather than two.
    """
    path = _write_fixture(tmp_path, rows=ROWS_DUPLICATES, unique_unread=True)
    outdir = tmp_path / "plots"

    run_evaluation_controls_from_path(str(path), str(outdir))

    assert "direct_target_rows=3" in (outdir / SUMMARY).read_text()


def test_unscorable_rows_are_reported_and_counted(tmp_path, capsys):
    """The dropped-row accounting the QC summary reads is unchanged."""
    assert perform_binary_evaluation(
        [1, 1, 0, 0], [0.01, np.nan, 0.4, 0.5], tmp_path, plot=False
    )

    assert "Dropped 1 row(s) with a missing label or p-value" in capsys.readouterr().out
    summary = (tmp_path / SUMMARY).read_text()
    assert "direct_target_rows=1" in summary
    assert "non_targeting_rows=2" in summary
    assert "unscorable_rows_dropped=1" in summary


@pytest.mark.parametrize(
    "mapping",
    [
        {"a": "A", "b": "B", "c": "C"},   # one-to-one
        {"a": "A", "b": "A", "c": "C"},   # two guides, one target
        {"a": "A"},                       # guides absent from the mapping
        {},                               # nothing maps
    ],
)
def test_mapping_a_categorical_matches_series_map(mapping):
    """The per-category lookup has to agree with the per-row one it replaced."""
    values = pd.Categorical(["a", "b", None, "c", "a", "b"])

    mapped = ec._map_categorical(values, mapping)

    expected = pd.Series(values).map(mapping)
    pd.testing.assert_series_equal(
        pd.Series(mapped).astype(object), expected.astype(object)
    )


def test_cli_writes_the_plots_the_task_publishes(tmp_path):
    path = _write_fixture(tmp_path, rows=ROWS)
    outdir = tmp_path / "plots"

    subprocess.run(
        [sys.executable, str(BIN_DIR / "evaluate_controls.py"), str(path),
         "--outdir", str(outdir)],
        check=True,
        capture_output=True,
        env={"PATH": "/usr/bin:/bin", "MPLBACKEND": "Agg",
             "MPLCONFIGDIR": str(tmp_path / "mpl"),
             "PYTHONPATH": str(BIN_DIR)},
    )

    assert sorted(p.name for p in outdir.iterdir()) == sorted(
        [SUMMARY, *PLOTS, *COMPOSITION]
    )
