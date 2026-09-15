"""``igv.py`` writes pipeline artifacts, so the fast path has to write the same bytes.

``evaluation_plot`` walked ``global_analysis_per_element_results`` -- 13,140,543
rows on the TAP-seq chr8 screen -- with ``DataFrame.iterrows`` and took 1.19 h to
emit 78,690 rows of bedpe/bedgraph. The rewrite is column-at-a-time, and these
tests pin the two things that are easy to lose on the way:

* ``to_csv`` formats from dtype. Gene coordinates arrive as float64 and guide
  coordinates as nullable Int64, so the same output column prints ``1000.0`` or
  ``11000`` depending on which dictionary entry fed it. A rewrite that routed
  coordinates through ``Series.map`` would unify those and silently rewrite
  every start/end in the file.
* ``row[a] == row[b]`` is Python ``==``: ``None == None`` is True. pandas'
  ``Series == Series`` forces False whenever either side is null, so the naive
  vectorisation drops rows the row loop kept.

Both the golden bytes below and the comparison against the pre-change script
cover those. GTF parsing is stubbed rather than run through ``gtfparse``: the
subject is the row loop, and a stub keeps the golden bytes independent of which
parser version happens to be installed.
"""

import gzip
import importlib.util
import os
import pathlib
import subprocess
import sys
import types

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

GLOBAL_KEY = "global_analysis_per_element_results"
LOCAL_KEY = "local_analysis_per_element_results"


# --------------------------------------------------------------------------- #
# A stand-in for gtfparse.read_gtf
# --------------------------------------------------------------------------- #

STUB_SOURCE = '''
"""Minimal stand-in for gtfparse, for the igv.py tests."""

import gzip
import re

import pandas as pd


class _Frame:
    def __init__(self, frame):
        self._frame = frame

    def to_pandas(self):
        return self._frame


def read_gtf(path, *args, **kwargs):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9:
                continue
            attributes = dict(re.findall(r'(\\S+) "([^"]*)"', fields[8]))
            rows.append(
                {
                    "seqname": fields[0],
                    "feature": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "gene_id": attributes.get("gene_id", ""),
                    # An absent gene_name stays absent: it is how a GTF row
                    # produces a null gtf_gene_name without going through the
                    # merge, which is the only way to reach `None == None`.
                    "gene_name": attributes.get("gene_name"),
                }
            )
    return _Frame(pd.DataFrame(rows))
'''

_stub_namespace = {}
exec(compile(STUB_SOURCE, "<gtfparse stub>", "exec"), _stub_namespace)
stub_read_gtf = _stub_namespace["read_gtf"]


def _import_igv(name, source):
    """Import an ``igv.py`` source with ``gtfparse`` stubbed out.

    ``igv.py`` does ``from gtfparse import read_gtf`` at import time, so the stub
    has to be in ``sys.modules`` for the import and is taken back out afterwards
    -- the module keeps its own reference, and nothing else in the suite should
    inherit a fake gtfparse.
    """
    path = pathlib.Path(source)
    stub = types.ModuleType("gtfparse")
    stub.read_gtf = stub_read_gtf
    previous = sys.modules.get("gtfparse")
    sys.modules["gtfparse"] = stub
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            sys.modules.pop("gtfparse", None)
        else:
            sys.modules["gtfparse"] = previous
    module.read_gtf = stub_read_gtf
    return module


igv = _import_igv("igv_current", BIN_DIR / "igv.py")


# --------------------------------------------------------------------------- #
# The pre-change script, found in history rather than copied
# --------------------------------------------------------------------------- #

ROW_LOOP_MARKER = "for index, row in merged_df.iterrows():"


def _git(*args):
    return subprocess.run(
        ["git", "-C", str(REPO_ROOT), *args],
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def _original_source():
    """The newest ``bin/igv.py`` in history that still used the row loop."""
    try:
        revisions = _git("rev-list", "HEAD", "--", "bin/igv.py").split()
    except (OSError, subprocess.CalledProcessError):
        return None
    for revision in revisions:
        try:
            source = _git("show", f"{revision}:bin/igv.py")
        except subprocess.CalledProcessError:
            continue
        if ROW_LOOP_MARKER in source:
            return source
    return None


ORIGINAL_SOURCE = _original_source()
needs_original = pytest.mark.skipif(
    ORIGINAL_SOURCE is None,
    reason="no row-loop bin/igv.py reachable in git history to compare against",
)


@pytest.fixture(scope="module")
def original(tmp_path_factory):
    if ORIGINAL_SOURCE is None:
        pytest.skip("no row-loop bin/igv.py in history")
    path = tmp_path_factory.mktemp("original") / "igv_original.py"
    path.write_text(ORIGINAL_SOURCE)
    return _import_igv("igv_original", path)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #

# gene_id -> gene_name, in file order. ENSG002 appears twice so
# drop_duplicates(subset=['gene_id2'], keep='first') has something to keep.
GTF_LINES = [
    ("chr8", "ENSG001.1", "SYMB1"),
    ("chr8", "ENSG002.5", "SYMB2"),
    ("chr9", "ENSG003.1", "SYMB3"),
    ("chr8", "ENSG002.9", "SYMB2_ALT"),
    ("chrX", "ENSG006.1", "ELEM2"),
    # No gene_name attribute at all: gtf_gene_name comes back null.
    ("chr8", "ENSG007.1", None),
]


def _write_gtf(path, lines=GTF_LINES, compress=False):
    text = ""
    for chrom, gene_id, gene_name in lines:
        attributes = f'gene_id "{gene_id}";'
        if gene_name is not None:
            attributes += f' gene_name "{gene_name}";'
        text += f"{chrom}\tTEST\tgene\t1\t2\t.\t+\t.\t{attributes}\n"
    if compress:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        pathlib.Path(path).write_text(text)
    return str(path)


def _anndata(var, n_obs=4):
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(n_obs)])
    return ad.AnnData(
        X=sparse.csr_matrix((n_obs, len(var)), dtype=np.float32),
        obs=obs,
        var=var,
    )


def _gene_var():
    """Gene coordinates as the screen stores them: float64, categorical chr."""
    return pd.DataFrame(
        {
            "gene_chr": pd.Categorical(["8", "8", "9", "9", "8"]),
            "gene_start": [1000.0, 3000.0, np.nan, 7000.0, 9000.0],
            "gene_end": [2000.0, 4000.0, 6000.0, np.nan, 9500.0],
        },
        index=pd.Index(
            ["ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005"], name="gene_id"
        ),
    )


def _guide_var():
    """Guide coordinates as the screen stores them: nullable Int64."""
    return pd.DataFrame(
        {
            "intended_target_name": pd.Categorical(
                [
                    "ELEM1",
                    "ELEM1",  # duplicate target: the first row's coords win
                    "non-targeting",  # skipped outright
                    "SYMB1",
                    "ENSG002",  # collides with a gene entry, which keeps priority
                    "ELEM2",
                ]
            ),
            "intended_target_chr": pd.Categorical(
                ["chr8", "chr8", "chrX", "chr8", "chr8", "chrX"]
            ),
            "intended_target_start": pd.array(
                [11000, 99999, 1, 13000, 15000, 17000], dtype="Int64"
            ),
            "intended_target_end": pd.array(
                [12000, 99999, 2, 14000, 16000, 18000], dtype="Int64"
            ),
        },
        index=pd.Index([f"g{i}" for i in range(6)], name="guide_id"),
    )


# gene_id, intended_target_name, log2_fc, p_value, and what each row exercises.
GLOBAL_ROWS = [
    ("ENSG001", "SYMB1", 0.5, 0.01),  # target == gtf name, placed -> bedgraph
    ("ENSG001", "ELEM1", -1.5, 0.002),  # differ, both placed -> bedpe
    ("ENSG002", "ELEM1", 0.25, 0.3),  # bedpe
    ("ENSG002", "ELEM2", np.nan, 0.4),  # dropna on log2_fc
    ("ENSG002", "ELEM2", 0.7, np.nan),  # dropna on p_value
    ("ENSG003", "ELEM1", 1.0, 0.05),  # gene has no coords -> dropped
    ("ENSG005", "ELEM2", -0.2, 0.6),  # gene not in the GTF -> null name -> bedpe
    ("ENSG001", "UNPLACED", 2.0, 0.07),  # target has no coords -> dropped
    ("ENSG002", "SYMB2", 0.9, 0.08),  # names match but target unplaced -> dropped
    ("ENSG006", "ELEM2", 0.33, 0.09),  # target == gtf name via ENSG006 -> bedgraph
    ("ENSG002", "SYMB2_ALT", 1.1, 0.11),  # the dropped duplicate GTF name
    ("ENSG002", "ELEM2", 0.44, 0.12),  # bedpe
]

LOCAL_ROWS = [
    ("ENSG001", "ELEM1", -0.75, 0.001),
    ("ENSG001", "SYMB1", 0.6, 0.02),
]


def _results_frame(
    rows, prefixes=("perturbo",), extra_columns=True, log2_fc_dtype=np.float64
):
    gene_id, target, log2_fc, p_value = (list(column) for column in zip(*rows))
    frame = {
        "gene_id": pd.Categorical(gene_id),
        "intended_target_name": pd.Categorical(target),
    }
    for prefix in prefixes:
        frame[f"{prefix}_log2_fc"] = np.array(log2_fc, dtype=log2_fc_dtype)
        frame[f"{prefix}_p_value"] = np.array(p_value, dtype=np.float64)
        if extra_columns:
            # Columns igv.py never reads; the real table is twenty-odd wide.
            frame[f"{prefix}_q_value"] = np.linspace(0, 1, len(rows))
            frame[f"{prefix}_negLog10p"] = np.linspace(1, 2, len(rows))
    if extra_columns:
        frame["intended_target_chr"] = pd.Categorical(["chr8"] * len(rows))
        frame["nPerturbedCells"] = pd.array(range(len(rows)), dtype="Int64")
    return pd.DataFrame(frame)


def _screen_mudata():
    """The dtype layout of the real screen: float64 genes, Int64 guides."""
    mdata = mu.MuData(
        {"gene": _anndata(_gene_var()), "guide": _anndata(_guide_var())}
    )
    mdata.uns[LOCAL_KEY] = _results_frame(LOCAL_ROWS, ("sceptre", "perturbo"))
    mdata.uns[GLOBAL_KEY] = _results_frame(GLOBAL_ROWS, ("perturbo",))
    mdata.uns["test_results"] = pd.DataFrame(
        {
            "gene_id": pd.Categorical([row[0] for row in GLOBAL_ROWS]),
            "intended_target_name": pd.Categorical([row[1] for row in GLOBAL_ROWS]),
            "log2_fc": np.array([row[2] for row in GLOBAL_ROWS], dtype=np.float64),
            "p_value": np.array([row[3] for row in GLOBAL_ROWS], dtype=np.float64),
        }
    )
    return mdata


def _edge_mudata():
    """Nulls, missing guide coordinates, and a target keyed to a gene entry.

    ``intended_target_name`` is object rather than categorical here so a literal
    ``None`` survives: paired with a GTF row that has no ``gene_name``, it is the
    ``None == None`` case where Python ``==`` and pandas ``==`` disagree.
    """
    gene_var = pd.DataFrame(
        {
            "gene_chr": pd.Categorical(["8", "8"]),
            "gene_start": [1000.0, 3000.0],
            "gene_end": [2000.0, 4000.0],
        },
        index=pd.Index(["ENSG001", "ENSG007"], name="gene_id"),
    )
    guide_var = pd.DataFrame(
        {
            # np.array, not pd.Series: a Series carries its own RangeIndex and
            # DataFrame(..., index=...) would reindex it to all-NaN.
            "intended_target_name": np.array(
                ["ELEM1", None, "ELEM_NA", "ENSG001"], dtype=object
            ),
            "intended_target_chr": np.array(["chr8"] * 4, dtype=object),
            # A missing guide coordinate: the guide branch has no nan guard, so
            # the pd.NA reaches the output column and the dtype turns object.
            "intended_target_start": pd.array(
                [11000, 21000, None, 31000], dtype="Int64"
            ),
            "intended_target_end": pd.array(
                [12000, 22000, None, 32000], dtype="Int64"
            ),
        },
        index=pd.Index([f"g{i}" for i in range(4)], name="guide_id"),
    )
    rows = [
        # ENSG007 has no gene_name in the GTF -> null name; the target is None
        # too, so Python `==` says these match and the row is a promoter.
        ("ENSG007", None, 0.5, 0.01),
        # np.nan against the same null name: `nan == nan` is False, so this one
        # takes the enhancer branch instead.
        ("ENSG007", np.nan, 0.6, 0.02),
        ("ENSG001", "ELEM1", -1.0, 0.03),
        # Target keyed to a gene entry: float64 coordinates land in the same
        # start1 column as ELEM1's Int64 ones.
        ("ENSG007", "ENSG001", 0.7, 0.04),
        # Target with missing coordinates: pd.NA in start1/end1.
        ("ENSG007", "ELEM_NA", 0.8, 0.05),
    ]
    gene_ids, targets, log2_fc, p_values = (list(column) for column in zip(*rows))
    mdata = mu.MuData({"gene": _anndata(gene_var), "guide": _anndata(guide_var)})
    mdata.uns[GLOBAL_KEY] = pd.DataFrame(
        {
            "gene_id": np.array(gene_ids, dtype=object),
            "intended_target_name": np.array(targets, dtype=object),
            "perturbo_log2_fc": np.array(log2_fc, dtype=np.float64),
            "perturbo_p_value": np.array(p_values, dtype=np.float64),
        }
    )
    return mdata


def _nothing_survives_mudata():
    """Every row filtered out, so both files must come out empty."""
    mdata = mu.MuData(
        {"gene": _anndata(_gene_var()), "guide": _anndata(_guide_var())}
    )
    mdata.uns[GLOBAL_KEY] = _results_frame(
        [
            ("ENSG003", "UNPLACED", 1.0, 0.1),  # neither side placed
            ("ENSG001", "ELEM1", np.nan, np.nan),  # dropna
        ],
        ("perturbo",),
    )
    return mdata


def _float32_metric_mudata():
    """What production actually stores: a float32 log2_fc.

    ``global_analysis_per_element_results`` keeps ``perturbo_log2_fc`` as
    float32 while ``perturbo_p_value`` stays float64. The row loop read those
    cells out of an object array, so each float32 arrived as a Python float and
    the output column inferred back to float64 -- the file carries the full
    float64 repr of the float32 value. Handing the float32 array straight to
    ``pd.DataFrame`` instead keeps it float32 and truncates every log2_fc in the
    file, which is how this went wrong the first time.
    """
    mdata = _screen_mudata()
    mdata.uns[GLOBAL_KEY] = _results_frame(
        GLOBAL_ROWS, ("perturbo",), log2_fc_dtype=np.float32
    )
    return mdata


MUDATA_BUILDERS = {
    "screen": _screen_mudata,
    "float32_metric": _float32_metric_mudata,
    "edge": _edge_mudata,
    "nothing_survives": _nothing_survives_mudata,
}


def _as_csv(frame):
    return frame.to_csv(sep="\t", index=False, header=False)


# --------------------------------------------------------------------------- #
# Golden bytes: what the row loop produced, spelled out
# --------------------------------------------------------------------------- #

# start1/end1 come from Int64 guide coordinates and print bare; start2/end2 come
# from float64 gene coordinates and print with a trailing .0. That asymmetry is
# the file format, not an accident.
EXPECTED_GLOBAL_BEDPE = (
    "chr8\t11000\t12000\t8\t1000.0\t2000.0\t0.002\t-1.5\n"
    "chr8\t11000\t12000\t8\t3000.0\t4000.0\t0.3\t0.25\n"
    "chrX\t17000\t18000\t8\t9000.0\t9500.0\t0.6\t-0.2\n"
    "chrX\t17000\t18000\t8\t3000.0\t4000.0\t0.12\t0.44\n"
)

EXPECTED_GLOBAL_BEDGRAPH = (
    "chr8\t13000\t14000\t0.01\t0.5\n" "chrX\t17000\t18000\t0.09\t0.33\n"
)

EXPECTED_LOCAL_BEDPE = "chr8\t11000\t12000\t8\t1000.0\t2000.0\t0.001\t-0.75\n"
EXPECTED_LOCAL_BEDGRAPH = "chr8\t13000\t14000\t0.02\t0.6\n"


def test_screen_bedpe_and_bedgraph_bytes(tmp_path):
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(_screen_mudata(), gtf, "perturbo", GLOBAL_KEY)
    assert _as_csv(bedpe) == EXPECTED_GLOBAL_BEDPE
    assert _as_csv(bedgraph) == EXPECTED_GLOBAL_BEDGRAPH


def test_coordinate_dtypes_are_not_unified(tmp_path):
    """The whole point of the golden bytes, asserted as dtypes."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(_screen_mudata(), gtf, "perturbo", GLOBAL_KEY)
    assert bedpe["start1"].dtype == np.dtype("int64")
    assert bedpe["end1"].dtype == np.dtype("int64")
    assert bedpe["start2"].dtype == np.dtype("float64")
    assert bedpe["end2"].dtype == np.dtype("float64")
    assert bedgraph["start"].dtype == np.dtype("int64")
    assert bedgraph["end"].dtype == np.dtype("int64")


def test_float32_metric_widens_to_float64(tmp_path):
    """A float32 log2_fc must print at float64 width, as the row loop made it.

    This is the one difference real data caught: the global table's log2_fc is
    float32, and keeping it float32 rewrote all 70,776 log2_fc values in the
    screen's bedpe (``-0.17034076`` for ``-0.17034076154232025``).
    """
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(
        _float32_metric_mudata(), gtf, "perturbo", GLOBAL_KEY
    )
    assert bedpe["log2_fc"].dtype == np.dtype("float64")
    assert bedgraph["log2_fc"].dtype == np.dtype("float64")
    assert _as_csv(bedpe) == (
        "chr8\t11000\t12000\t8\t1000.0\t2000.0\t0.002\t-1.5\n"
        "chr8\t11000\t12000\t8\t3000.0\t4000.0\t0.3\t0.25\n"
        "chrX\t17000\t18000\t8\t9000.0\t9500.0\t0.6\t-0.20000000298023224\n"
        "chrX\t17000\t18000\t8\t3000.0\t4000.0\t0.12\t0.4399999976158142\n"
    )
    assert _as_csv(bedgraph) == (
        "chr8\t13000\t14000\t0.01\t0.5\n"
        "chrX\t17000\t18000\t0.09\t0.33000001311302185\n"
    )


def test_column_order(tmp_path):
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(_screen_mudata(), gtf, "perturbo", GLOBAL_KEY)
    assert list(bedpe.columns) == [
        "chr1",
        "start1",
        "end1",
        "chr2",
        "start2",
        "end2",
        "p_value",
        "log2_fc",
    ]
    assert list(bedgraph.columns) == ["chr", "start", "end", "p_value", "log2_fc"]


def test_generic_columns_path(tmp_path):
    """``method=None`` reads log2_fc/p_value and must agree with the method path."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(_screen_mudata(), gtf, None, "test_results")
    assert _as_csv(bedpe) == EXPECTED_GLOBAL_BEDPE
    assert _as_csv(bedgraph) == EXPECTED_GLOBAL_BEDGRAPH


def test_null_name_pair_takes_the_promoter_branch(tmp_path):
    """``None == None`` is True, which pandas' ``Series == Series`` denies."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(_edge_mudata(), gtf, "perturbo", GLOBAL_KEY)
    # Row 0's None/None pair matches, so it is a promoter at the coordinates the
    # guide row with a None target contributed. A pandas-`==` rewrite would send
    # it to the enhancer branch and this file would come out empty.
    assert _as_csv(bedgraph) == "chr8\t21000\t22000\t0.01\t0.5\n"
    assert len(bedgraph) == 1
    # Row 1's nan does not match the same null name, and nan is not the key the
    # None-targeted guide created, so it reaches neither file: 3 of 5 rows left.
    assert len(bedpe) == 3


def test_missing_guide_coordinates_reach_the_output(tmp_path):
    """The guide branch has no nan guard, so pd.NA is part of the contract."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, _ = igv.igv(_edge_mudata(), gtf, "perturbo", GLOBAL_KEY)
    # Two coordinate sources feed start1 here -- Int64 guide entries, float64
    # gene entries and a pd.NA -- so the column is object and the missing pair
    # writes as empty fields rather than being dropped or filled.
    assert bedpe["start1"].dtype == np.dtype("object")
    assert _as_csv(bedpe).splitlines()[-1] == "chr8\t\t\t8\t3000.0\t4000.0\t0.05\t0.8"
    # The gene-keyed target keeps its float64 coordinates in the same column.
    assert _as_csv(bedpe).splitlines()[1] == "8\t1000.0\t2000.0\t8\t3000.0\t4000.0\t0.04\t0.7"


def test_empty_results_write_empty_files(tmp_path):
    gtf = _write_gtf(tmp_path / "genes.gtf")
    bedpe, bedgraph = igv.igv(
        _nothing_survives_mudata(), gtf, "perturbo", GLOBAL_KEY
    )
    assert bedpe.empty and bedgraph.empty
    assert _as_csv(bedpe) == ""
    assert _as_csv(bedgraph) == ""


def test_duplicate_gtf_gene_id_keeps_the_first_name(tmp_path):
    """ENSG002 is in the GTF twice; only SYMB2 may be used, never SYMB2_ALT."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    mdata = _screen_mudata()
    mdata.uns[GLOBAL_KEY] = _results_frame(
        [
            ("ENSG002", "SYMB2", 1.0, 0.1),  # matches the kept name
            ("ENSG002", "SYMB2_ALT", 2.0, 0.2),  # matches the dropped one
        ],
        ("perturbo",),
    )
    # Give both names coordinates so only the name resolution decides.
    guide_var = mdata.mod["guide"].var.copy()
    guide_var["intended_target_name"] = pd.Categorical(
        ["SYMB2", "SYMB2_ALT", "non-targeting", "SYMB1", "ENSG002", "ELEM2"]
    )
    mdata.mod["guide"].var = guide_var
    bedpe, bedgraph = igv.igv(mdata, gtf, "perturbo", GLOBAL_KEY)
    # SYMB2 is the surviving GTF name, so row 0 is a promoter and row 1 is not.
    assert len(bedgraph) == 1
    assert len(bedpe) == 1
    assert _as_csv(bedgraph) == "chr8\t11000\t12000\t0.1\t1.0\n"


def test_process_coordinates_priority_and_nan_guard():
    coordinate_dict = igv.process_coordinates(_screen_mudata())
    assert set(coordinate_dict) == {
        "ENSG001",
        "ENSG002",
        "ENSG005",
        "ELEM1",
        "SYMB1",
        "ELEM2",
    }
    # Genes with a nan start or end are skipped; "non-targeting" never lands.
    assert "ENSG003" not in coordinate_dict
    assert "ENSG004" not in coordinate_dict
    assert "non-targeting" not in coordinate_dict
    # The gene entry wins over the guide row that reuses its name, and the
    # first guide row for a repeated target wins over later ones.
    assert coordinate_dict["ENSG002"] == ["8", 3000.0, 4000.0]
    assert coordinate_dict["ELEM1"][1] == 11000
    # Scalar types, not just values: they decide the output dtype, and these are
    # plain Python scalars because iterrows reads DataFrame.values, which casts a
    # mixed .var frame to object and unboxes every numpy scalar in it.
    assert type(coordinate_dict["ENSG002"][1]) is float
    assert type(coordinate_dict["ELEM1"][1]) is int


# --------------------------------------------------------------------------- #
# Against the pre-change script
# --------------------------------------------------------------------------- #


@needs_original
@pytest.mark.parametrize("builder", sorted(MUDATA_BUILDERS))
@pytest.mark.parametrize("method", ["perturbo", None])
def test_matches_row_loop_output(original, builder, method, tmp_path):
    """Same bytes and same dtypes as the row loop, table for table."""
    gtf = _write_gtf(tmp_path / "genes.gtf")
    key = GLOBAL_KEY
    if method is None:
        mdata_old, mdata_new = MUDATA_BUILDERS[builder](), MUDATA_BUILDERS[builder]()
        for mdata in (mdata_old, mdata_new):
            table = mdata.uns[key].rename(
                columns={
                    "perturbo_log2_fc": "log2_fc",
                    "perturbo_p_value": "p_value",
                }
            )
            mdata.uns[key] = table
    else:
        mdata_old, mdata_new = MUDATA_BUILDERS[builder](), MUDATA_BUILDERS[builder]()

    old_bedpe, old_bedgraph = original.igv(mdata_old, gtf, method, key)
    new_bedpe, new_bedgraph = igv.igv(mdata_new, gtf, method, key)

    assert _as_csv(new_bedpe) == _as_csv(old_bedpe)
    assert _as_csv(new_bedgraph) == _as_csv(old_bedgraph)
    assert list(new_bedpe.columns) == list(old_bedpe.columns)
    assert list(new_bedgraph.columns) == list(old_bedgraph.columns)
    assert new_bedpe.dtypes.to_dict() == old_bedpe.dtypes.to_dict()
    assert new_bedgraph.dtypes.to_dict() == old_bedgraph.dtypes.to_dict()


@needs_original
@pytest.mark.parametrize("builder", sorted(MUDATA_BUILDERS))
def test_process_coordinates_matches_row_loop(original, builder):
    old = original.process_coordinates(MUDATA_BUILDERS[builder]())
    new = igv.process_coordinates(MUDATA_BUILDERS[builder]())
    assert list(new) == list(old)  # same keys, same insertion order
    for key in old:
        for mine, theirs in zip(new[key], old[key]):
            assert type(mine) is type(theirs)
            assert mine is theirs or mine == theirs or (
                pd.isna(mine) and pd.isna(theirs)
            )


# --------------------------------------------------------------------------- #
# End to end, as the module invokes it
# --------------------------------------------------------------------------- #


def _run_script(script, mudata_path, gtf, workdir, extra_args):
    """Run igv.py the way modules/local/evaluation_plot does."""
    stub_dir = workdir.parent / "gtfparse_stub"
    stub_dir.mkdir(exist_ok=True)
    (stub_dir / "gtfparse.py").write_text(STUB_SOURCE)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(stub_dir), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    workdir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [sys.executable, str(script), str(mudata_path), "--gtf", gtf, *extra_args],
        cwd=str(workdir),
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    output = workdir / "evaluation_output"
    return {
        path.name: path.read_bytes() for path in sorted(output.iterdir())
    }, result.stdout


@needs_original
@pytest.mark.parametrize(
    "extra_args",
    [["--default"], ["--results_key", "test_results"]],
    ids=["default", "single_key"],
)
def test_script_writes_identical_files(original, tmp_path, extra_args):
    mudata_path = tmp_path / "inference_mudata.h5mu"
    _screen_mudata().write(mudata_path)
    # .gz, as the module stages it.
    gtf = _write_gtf(tmp_path / "gencode_gtf.gtf.gz", compress=True)

    original_script = tmp_path / "igv_original_script.py"
    original_script.write_text(ORIGINAL_SOURCE)

    old_files, _ = _run_script(
        original_script, mudata_path, gtf, tmp_path / "old", extra_args
    )
    new_files, _ = _run_script(
        BIN_DIR / "igv.py", mudata_path, gtf, tmp_path / "new", extra_args
    )

    assert sorted(new_files) == sorted(old_files)
    assert new_files == old_files, "evaluation_output bytes differ"
    for name, payload in new_files.items():
        assert payload, f"{name} unexpectedly empty"


def _missing_category_mudata():
    """A result table with a missing target name and a missing gene_id.

    Written to h5mu these are categorical code -1, which the slim reader decodes
    as ``None`` and ``read_h5mu`` decodes as nan -- and Python ``==`` tells those
    two apart. igv.py normalises, so the branch does not depend on the reader.
    """
    mdata = _screen_mudata()
    frame = _results_frame(GLOBAL_ROWS, ("perturbo",))
    target = frame["intended_target_name"].astype(object).to_numpy()
    gene_id = frame["gene_id"].astype(object).to_numpy()
    target[0] = None
    gene_id[2] = None
    frame["intended_target_name"] = pd.Categorical(target)
    frame["gene_id"] = pd.Categorical(gene_id)
    mdata.uns[GLOBAL_KEY] = frame
    return mdata


def test_read_results_normalizes_the_slim_readers_nulls(tmp_path):
    mudata_path = tmp_path / "inference_mudata.h5mu"
    _missing_category_mudata().write(mudata_path)
    frame = igv._read_results(
        str(mudata_path),
        GLOBAL_KEY,
        ["gene_id", "intended_target_name", "perturbo_log2_fc", "perturbo_p_value"],
    )
    # Only what was asked for, and no None left where read_h5mu would give nan.
    assert list(frame.columns) == [
        "gene_id",
        "intended_target_name",
        "perturbo_log2_fc",
        "perturbo_p_value",
    ]
    assert frame["intended_target_name"].isna().sum() == 1
    assert frame["gene_id"].isna().sum() == 1
    assert not any(value is None for value in frame["intended_target_name"])
    assert not any(value is None for value in frame["gene_id"])
    # Same values read_h5mu gives for those columns.
    expected = mu.read_h5mu(mudata_path).uns[GLOBAL_KEY]
    for column in frame.columns:
        got = [None if pd.isna(v) else v for v in frame[column]]
        want = [None if pd.isna(v) else v for v in expected[column]]
        assert got == want, column


@needs_original
def test_script_identical_with_missing_categories(tmp_path):
    """The whole script, on a table the two readers decode differently."""
    mudata_path = tmp_path / "inference_mudata.h5mu"
    _missing_category_mudata().write(mudata_path)
    gtf = _write_gtf(tmp_path / "gencode_gtf.gtf.gz", compress=True)
    original_script = tmp_path / "igv_original_script.py"
    original_script.write_text(ORIGINAL_SOURCE)
    old_files, _ = _run_script(
        original_script, mudata_path, gtf, tmp_path / "old", ["--default"]
    )
    new_files, _ = _run_script(
        BIN_DIR / "igv.py", mudata_path, gtf, tmp_path / "new", ["--default"]
    )
    assert new_files == old_files


def test_default_writes_the_expected_file_set(tmp_path):
    mudata_path = tmp_path / "inference_mudata.h5mu"
    _screen_mudata().write(mudata_path)
    gtf = _write_gtf(tmp_path / "gencode_gtf.gtf.gz", compress=True)
    files, stdout = _run_script(
        BIN_DIR / "igv.py", mudata_path, gtf, tmp_path / "run", ["--default"]
    )
    assert sorted(files) == [
        "global_analysis_perturbo.bedgraph",
        "global_analysis_perturbo.bedpe",
        "local_analysis_perturbo.bedgraph",
        "local_analysis_perturbo.bedpe",
        "local_analysis_sceptre.bedgraph",
        "local_analysis_sceptre.bedpe",
    ]
    assert files["global_analysis_perturbo.bedpe"].decode() == EXPECTED_GLOBAL_BEDPE
    assert (
        files["global_analysis_perturbo.bedgraph"].decode()
        == EXPECTED_GLOBAL_BEDGRAPH
    )
    assert files["local_analysis_perturbo.bedpe"].decode() == EXPECTED_LOCAL_BEDPE
    assert (
        files["local_analysis_perturbo.bedgraph"].decode() == EXPECTED_LOCAL_BEDGRAPH
    )
    assert "Available methods" in stdout


def test_missing_results_key_is_a_warning_not_a_failure(tmp_path):
    mdata = mu.MuData(
        {"gene": _anndata(_gene_var()), "guide": _anndata(_guide_var())}
    )
    mudata_path = tmp_path / "inference_mudata.h5mu"
    mdata.write(mudata_path)
    gtf = _write_gtf(tmp_path / "genes.gtf")
    files, stdout = _run_script(
        BIN_DIR / "igv.py", mudata_path, gtf, tmp_path / "run", ["--default"]
    )
    assert files == {}
    assert "not found in mdata.uns" in stdout
