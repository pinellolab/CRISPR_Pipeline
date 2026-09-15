"""The categorical intended-target mask must agree with the row-wise one.

``intended_target.py`` used to map ``intended_target_name`` onto a copy of the
whole result table and let ``direct_target_mask`` normalize identifiers once per
row. On the TAP-seq chr8 screen that table has 52,006,760 rows drawn from 4,120
guides and 1,041 intended targets, so the normalization -- strip, upper, an
Ensembl-version regex, blank-to-missing -- ran thousands of times per distinct
value: 13.34 s at 2M rows, 52.17 s at 8M, 131.78 s at 20M.

It now normalizes per category instead. That is only worth having if the mask is
the same one, so the reference expressions below are the code that was replaced,
and the tests assert the two agree -- on every identifier shape the old path
handled specially, and on a table large enough that a rare category matters.
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

import intended_target as it
from inference_target_matching import direct_target_mask


# ---------------------------------------------------------------------
# The replaced implementation, kept verbatim as the reference
# ---------------------------------------------------------------------
def reference_mask(trans_results, intended_map):
    results = trans_results.copy()
    results["intended_target_name"] = results["guide_id"].map(intended_map)
    return direct_target_mask(results)


def reference_filter(
    trans_results,
    guide_var,
    log2fc_col="log2_fc",
    pvalue_col="p_value",
    non_targeting_name="non-targeting",
):
    guide_meta = it.build_guide_meta(guide_var, non_targeting_name=non_targeting_name)
    intended_map = guide_meta.set_index("guide_id")["intended_target_name"].to_dict()

    results = trans_results.copy()
    results["intended_target_name"] = results["guide_id"].map(intended_map)

    intended = results[direct_target_mask(results)].copy()
    intended = intended.drop_duplicates(subset=["guide_id", "gene_id"])
    return guide_meta.merge(
        intended[["guide_id", "gene_id", log2fc_col, pvalue_col]],
        on="guide_id",
        how="left",
    )


def reference_evaluation_table(
    trans_results,
    guide_var,
    pvalue_col="p_value",
    non_targeting_name="non-targeting",
    random_state=42,
):
    guide_meta = it.build_guide_meta(guide_var, non_targeting_name=non_targeting_name)
    intended_map = guide_meta.set_index("guide_id")["intended_target_name"].to_dict()
    targeting_map = guide_meta.set_index("guide_id")["targeting"].to_dict()

    results = trans_results.copy()
    results["intended_target_name"] = results["guide_id"].map(intended_map)
    results["targeting"] = results["guide_id"].map(targeting_map)

    positives = results[
        direct_target_mask(results) & results[pvalue_col].notna()
    ].copy()
    positives["direct_target"] = 1

    all_target_gene_ids = positives["gene_id"].dropna().unique()

    negatives = results[
        (results["targeting"] == False) &
        (results["gene_id"].isin(all_target_gene_ids)) &
        (results[pvalue_col].notna())
    ].copy()
    negatives["direct_target"] = 0

    n_pos = len(positives)
    n_neg = len(negatives)
    if n_neg > n_pos and n_pos > 0:
        negatives = negatives.sample(n=n_pos, random_state=random_state)
    elif n_neg == 0 or n_pos == 0:
        return pd.DataFrame()

    eval_table = pd.concat([positives, negatives], ignore_index=True)
    eval_table = eval_table[
        ["guide_id", "gene_id", pvalue_col, "direct_target"]
    ].copy()
    return eval_table.rename(columns={pvalue_col: "p_value"})


# ---------------------------------------------------------------------
# Fixtures: one guide panel exercising every identifier shape
# ---------------------------------------------------------------------
# Each guide below is the shape it is named for. The genes they are tested
# against are in GENES, so a guide's row set covers matching and non-matching
# genes alike.
GUIDE_VAR = pd.DataFrame(
    {
        "guide_id": [
            "g_symbol",        # target is a gene symbol, matches gene_name
            "g_ensembl",       # target is an Ensembl id, matches gene_id
            "g_versioned",     # target carries a version suffix the regex strips
            "g_padded_case",   # target needs strip() and upper() to match
            "g_blank",         # target is "", which must not match another blank
            "g_spaces",        # target is whitespace only, likewise
            "g_missing_name",  # target is NaN
            "g_element",       # target is a genomic interval, matches nothing
            "nt_1",            # non-targeting
            "nt_2",            # non-targeting
            "g_untested",      # in guide.var, absent from the result table
        ],
        "intended_target_name": [
            "FAM83A",
            "ENSG00000147689",
            "ENSG00000104312.7",
            "  ripk2  ",
            "",
            "   ",
            np.nan,
            "chr8:127735434-127742951",
            "non-targeting",
            "non-targeting",
            "MYC",
        ],
        "targeting": [
            True, True, True, True, True, True, True, True, False, False, True,
        ],
        "gene_name": [
            "FAM83A", "FAM83A", "RIPK2", "RIPK2", "", "", "", "ELEM", "", "", "MYC",
        ],
    }
)
GUIDE_VAR.index = GUIDE_VAR["guide_id"]

# gene_id / gene_name pairs the guides are tested against. The blank pair is
# what upstream formatting leaves behind for a gene with no symbol, and the
# versioned id is the other half of the version-suffix case.
GENES = [
    ("ENSG00000147689", "FAM83A"),
    ("ENSG00000104312", "RIPK2"),
    ("ENSG00000104312.3", "RIPK2"),
    ("ENSG00000136997", "MYC"),
    ("", ""),
    (np.nan, np.nan),
    ("ENSG00000075624", "ACTB"),
]

# Guides in the result table but not in guide.var: the map yields NaN for them.
UNKNOWN_GUIDES = ["g_not_in_var", "g_also_not_in_var"]


def _results_frame(dup_rows=True, with_gene_name=True, nan_pvalues=True):
    guides = [g for g in GUIDE_VAR["guide_id"] if g != "g_untested"]
    guides = guides + UNKNOWN_GUIDES

    rows = []
    for guide in guides:
        for gene_id, gene_name in GENES:
            rows.append((guide, gene_id, gene_name))
    if dup_rows:
        # A repeated guide/gene pair, which drop_duplicates must collapse, and a
        # repeat of a matching pair so the surviving row is the first one.
        rows += [
            ("g_symbol", "ENSG00000147689", "FAM83A"),
            ("g_symbol", "ENSG00000147689", "FAM83A"),
            ("g_ensembl", "ENSG00000147689", "FAM83A"),
        ]

    frame = pd.DataFrame(rows, columns=["guide_id", "gene_id", "gene_name"])
    rng = np.random.default_rng(11)
    frame["log2_fc"] = rng.normal(size=len(frame))
    frame["p_value"] = rng.random(size=len(frame))
    if nan_pvalues:
        frame.loc[frame.index[::9], "p_value"] = np.nan
    # A 'targeting' column of its own, which the mapped one must shadow.
    frame["targeting"] = True
    if not with_gene_name:
        frame = frame.drop(columns=["gene_name"])
    return frame


def _categorical(frame):
    out = frame.copy()
    for column in ("guide_id", "gene_id", "gene_name"):
        if column in out.columns:
            out[column] = out[column].astype("category")
    return out


def _intended_map(guide_var=GUIDE_VAR):
    guide_meta = it.build_guide_meta(guide_var)
    return guide_meta.set_index("guide_id")["intended_target_name"].to_dict()


# ---------------------------------------------------------------------
# Mask equivalence
# ---------------------------------------------------------------------
@pytest.mark.parametrize("as_categorical", [False, True])
@pytest.mark.parametrize("with_gene_name", [True, False])
def test_mask_matches_row_wise_reference(as_categorical, with_gene_name):
    frame = _results_frame(with_gene_name=with_gene_name)
    if as_categorical:
        frame = _categorical(frame)

    expected = reference_mask(frame, _intended_map())
    actual = it._intended_target_mask(frame, _intended_map())

    assert actual.tolist() == expected.tolist()
    assert list(actual.index) == list(frame.index)
    # The panel must actually exercise both outcomes, or agreement is vacuous.
    assert expected.any() and not expected.all()


def test_mask_matches_reference_on_a_few_hundred_thousand_rows():
    frame = _results_frame()
    reps = 300_000 // len(frame) + 1
    big = pd.concat([frame] * reps, ignore_index=True)
    # Shuffle so the rare categories are not contiguous, then re-index oddly:
    # the mask is aligned by index and must not assume a RangeIndex.
    big = big.sample(frac=1.0, random_state=7).reset_index(drop=True)
    big.index = big.index * 3 + 1
    assert len(big) >= 300_000

    expected = reference_mask(big, _intended_map())
    actual = it._intended_target_mask(big, _intended_map())

    assert actual.tolist() == expected.tolist()
    assert list(actual.index) == list(big.index)
    assert expected.sum() > 0


def test_mask_raises_the_same_error_when_gene_id_is_absent():
    frame = _results_frame().drop(columns=["gene_id"])
    with pytest.raises(KeyError, match="gene_id"):
        it._intended_target_mask(frame, _intended_map())


def test_mask_is_empty_for_an_empty_result_table():
    frame = _results_frame().iloc[:0]
    expected = reference_mask(frame, _intended_map())
    actual = it._intended_target_mask(frame, _intended_map())
    assert actual.tolist() == expected.tolist() == []


def test_mask_handles_a_guide_id_column_that_is_entirely_missing():
    # No guide categories at all, so every row lands on the trailing missing
    # slot and the mapped target column has no categories to point at.
    frame = _results_frame()
    frame["guide_id"] = np.nan
    expected = reference_mask(frame, _intended_map())
    actual = it._intended_target_mask(frame, _intended_map())
    assert actual.tolist() == expected.tolist()
    assert not expected.any()


def test_mask_handles_a_guide_panel_with_no_usable_targets():
    guide_var = GUIDE_VAR[GUIDE_VAR["guide_id"].isin(["g_blank", "nt_1"])]
    frame = _results_frame()
    expected = reference_mask(frame, _intended_map(guide_var))
    actual = it._intended_target_mask(frame, _intended_map(guide_var))
    assert actual.tolist() == expected.tolist()
    assert not expected.any()


# ---------------------------------------------------------------------
# The two callers
# ---------------------------------------------------------------------
@pytest.mark.parametrize("as_categorical", [False, True])
def test_filter_to_intended_targets_matches_reference(as_categorical):
    frame = _results_frame()
    if as_categorical:
        frame = _categorical(frame)

    expected = reference_filter(frame, GUIDE_VAR)
    actual = it.filter_to_intended_targets(frame, GUIDE_VAR)

    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    # Guides that matched a gene carry a result; the rest are left NaN by the
    # left merge, which is what the metrics then drop.
    assert actual["log2_fc"].notna().any()
    assert actual["log2_fc"].isna().any()


@pytest.mark.parametrize("as_categorical", [False, True])
def test_build_evaluation_table_matches_reference(as_categorical):
    frame = _results_frame()
    if as_categorical:
        frame = _categorical(frame)

    expected = reference_evaluation_table(frame, GUIDE_VAR)
    actual = it.build_evaluation_table(frame, GUIDE_VAR)

    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    assert set(actual["direct_target"]) == {0, 1}


def test_build_evaluation_table_ignores_the_tables_own_targeting_column():
    # The result table carries targeting=True for every row; the negatives come
    # from guide.var's non-targeting guides, so they must survive.
    frame = _results_frame()
    assert frame["targeting"].all()
    actual = it.build_evaluation_table(frame, GUIDE_VAR)
    negatives = actual[actual["direct_target"] == 0]
    assert set(negatives["guide_id"]) <= {"nt_1", "nt_2"}
    assert len(negatives) > 0


def test_build_evaluation_table_downsamples_negatives_identically():
    # Many more negatives than positives, so sample() runs; it must pick the
    # same rows as it did when the frame still carried all its columns.
    extra_nt = pd.DataFrame(
        {
            "guide_id": [f"nt_x{i:02d}" for i in range(40)],
            "intended_target_name": ["non-targeting"] * 40,
            "targeting": [False] * 40,
            "gene_name": [""] * 40,
        }
    )
    extra_nt.index = extra_nt["guide_id"]
    guide_var = pd.concat([GUIDE_VAR, extra_nt])

    frame = _results_frame(nan_pvalues=False)
    extra = pd.DataFrame(
        [
            (guide, gene_id, gene_name)
            for guide in extra_nt["guide_id"]
            for gene_id, gene_name in GENES
        ],
        columns=["guide_id", "gene_id", "gene_name"],
    )
    rng = np.random.default_rng(29)
    extra["log2_fc"] = rng.normal(size=len(extra))
    extra["p_value"] = rng.random(size=len(extra))
    extra["targeting"] = True
    frame = pd.concat([frame, extra], ignore_index=True)

    expected = reference_evaluation_table(frame, guide_var)
    actual = it.build_evaluation_table(frame, guide_var)
    n_pos = int((expected["direct_target"] == 1).sum())
    assert n_pos > 0
    assert int((expected["direct_target"] == 0).sum()) == n_pos
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)


def test_build_evaluation_table_is_empty_without_positives():
    guide_var = GUIDE_VAR[GUIDE_VAR["guide_id"].isin(["g_element", "nt_1"])]
    frame = _results_frame()
    assert it.build_evaluation_table(frame, guide_var).empty
    assert reference_evaluation_table(frame, guide_var).empty


# ---------------------------------------------------------------------
# End to end, through the written artifacts
# ---------------------------------------------------------------------
def _write_mudata(path, results, guide_var):
    # h5 cannot hold a NaN in a variable-length string dataset, so the columns
    # that carry one go out categorical -- which is how mergedResults writes
    # low-cardinality text anyway, and it is the encoding the slim reader
    # decodes back to object with None for the missing codes.
    results = results.copy()
    for column in ("gene_id", "gene_name", "guide_id"):
        results[column] = results[column].astype("category")
    guide_var = guide_var.copy()
    guide_var["intended_target_name"] = guide_var["intended_target_name"].astype(
        "category"
    )

    obs = pd.DataFrame(index=[f"cell{i}" for i in range(4)])
    gene_ids = [g for g, _ in GENES if isinstance(g, str) and g]
    gene = ad.AnnData(
        X=sparse.csr_matrix(
            np.arange(4 * len(gene_ids), dtype=np.float32).reshape(4, len(gene_ids))
        ),
        obs=obs.copy(),
        var=pd.DataFrame(index=gene_ids),
    )
    guide = ad.AnnData(
        X=sparse.csr_matrix(np.eye(4, len(guide_var), dtype=np.float32)),
        obs=obs.copy(),
        var=guide_var.copy(),
    )
    guide.layers["guide_assignment"] = guide.X.copy()
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["global_analysis_per_guide_results"] = results
    mdata.write(path)
    return path


def test_run_intended_target_qc_writes_the_reference_metrics(tmp_path):
    """The written artifacts must match what the row-wise path produced.

    The expected metrics are computed here from the reference expressions
    above, so this pins the artifacts to the replaced implementation rather
    than to numbers typed out by hand.
    """
    results = _results_frame()
    path = _write_mudata(tmp_path / "inference_mudata.h5mu", results, GUIDE_VAR)
    outdir = tmp_path / "out"

    it.run_intended_target_qc(
        input_path=str(path),
        outdir=str(outdir),
        results_key="auto",
        log2fc_col="auto",
        pvalue_col="auto",
    )

    metrics = pd.read_csv(outdir / "intended_target_metrics.tsv", sep="\t")
    written = pd.read_csv(outdir / "intended_target_results.tsv", sep="\t")

    # Rebuild the expectation from the reference expressions, reading the same
    # slim-loaded frame the script reads.
    loaded = it.read_result_columns(
        path,
        "global_analysis_per_guide_results",
        list(it.QC_RESULT_COLUMNS) + ["log2_fc", "p_value"],
    )
    guide_var = mu.read_h5mu(path).mod["guide"].var
    expected_intended = reference_filter(loaded, guide_var)
    expected_metrics = it.compute_knockdown_metrics(expected_intended)
    expected_eval = reference_evaluation_table(loaded, guide_var)
    auroc, auprc, _ = it.compute_auroc_auprc(expected_eval)

    assert metrics["n_guides_total"].iloc[0] == expected_metrics["n_guides_total"]
    assert metrics["n_guides_tested"].iloc[0] == expected_metrics["n_guides_tested"]
    assert (
        metrics["n_strong_knockdowns"].iloc[0]
        == expected_metrics["n_strong_knockdowns"]
    )
    assert metrics["n_significant"].iloc[0] == expected_metrics["n_significant"]
    assert metrics["median_log2fc"].iloc[0] == pytest.approx(
        expected_metrics["median_log2fc"]
    )
    assert metrics["auroc"].iloc[0] == pytest.approx(auroc)
    assert metrics["auprc"].iloc[0] == pytest.approx(auprc)
    assert metrics["n_eval_positives"].iloc[0] == int(
        (expected_eval["direct_target"] == 1).sum()
    )
    assert metrics["n_eval_negatives"].iloc[0] == int(
        (expected_eval["direct_target"] == 0).sum()
    )

    # Every guide in guide.var gets a row, including the one never tested; a
    # guide that matched two genes (a versioned id and its unversioned form)
    # gets one row per match, as the left merge has always produced.
    assert set(written["guide_id"]) == set(GUIDE_VAR["guide_id"])
    assert len(written) == len(expected_intended)
    assert written.loc[written["guide_id"] == "g_untested", "log2_fc"].isna().all()

    for name in (
        "intended_target_volcano.png",
        "intended_target_log2fc_distribution.png",
        "intended_target_roc_pr_curves.png",
    ):
        assert (outdir / name).exists()
