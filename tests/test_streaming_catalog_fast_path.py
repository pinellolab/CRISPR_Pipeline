"""The enriched-Parquet catalog fast path must actually be taken.

Every failure inside ``try_write_enriched_parquet_catalog`` is swallowed and
demoted to the pandas path, which is correct but materialises the whole
guide-by-gene table -- tens of millions of rows on a real screen. A regression
there is therefore invisible except as runtime, so these tests assert on the
boolean the fast path returns rather than only on the table it writes.

The case that motivated them: ``perturbo_cis_*`` are the local table's
``perturbo_*`` columns under another name, produced by a rename in the pandas
path, so they exist under that name in neither input file. Adding them to
OUTPUT_COLUMNS without teaching the fast path the same rename made every
catalog fall back.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import build_catalog_per_element_output as element_builder
import build_catalog_per_guide_output as guide_builder
import streaming_catalog_io
from analysis_output_formatting import neg_log10

pytest.importorskip("polars")


def _streaming_sink_available():
    """Whether this Polars can sink a sorted plan at all.

    Polars 0.20 refuses ``sink_parquet(maintain_order=True)`` with "not yet
    supported in standard engine", so the fast path cannot be taken there
    regardless of the column wiring these tests are about. Skip rather than fail
    when the suite runs inside such an image.
    """
    import tempfile

    import polars as pl
    from polars_compat import SINK_ENGINE, SINK_ORDER

    with tempfile.TemporaryDirectory() as directory:
        source = pathlib.Path(directory) / "in.parquet"
        pl.DataFrame({"a": [2, 1]}).write_parquet(source)
        try:
            pl.scan_parquet(source).sort("a", **SINK_ORDER).sink_parquet(
                pathlib.Path(directory) / "out.parquet", **SINK_ENGINE
            )
        except Exception:
            return False
    return True


pytestmark = pytest.mark.skipif(
    not _streaming_sink_available(),
    reason="this Polars cannot sink a sorted plan; the fast path is unavailable",
)

PVALUE_FLOOR = 1e-300
# What mergedResults writes into the local table. Stated outright rather than
# read back from the builder's alias map, so these tests still describe the real
# input shape on a checkout whose fast path ignores these columns entirely.
LOCAL_PERTURBO_COLUMNS = [
    "perturbo_log2_fc",
    "perturbo_p_value",
    "perturbo_q_value",
    "perturbo_fc_se",
    "perturbo_negLog10p",
]
BUILDERS = [
    pytest.param(guide_builder, "build_catalog_per_guide_output", id="per-guide"),
    pytest.param(element_builder, "build_catalog_per_element_output", id="per-element"),
]


def _fast_path_kwargs(module, builder_attr):
    """The arguments the builder itself passes, captured rather than restated."""
    captured = {}
    original = module.try_write_enriched_parquet_catalog

    def spy(**kwargs):
        captured.update(kwargs)
        return True  # claim success so the builder stops before reading MuData

    module.try_write_enriched_parquet_catalog = spy
    try:
        getattr(module, builder_attr)("local", "global", "mudata", "out")
    finally:
        module.try_write_enriched_parquet_catalog = original
    return captured


def _frame(columns, n, pvalues):
    rng = np.random.default_rng(7)
    data = {}
    for column in columns:
        if column == "gene_id":
            data[column] = [f"ENSG{i:05d}" for i in range(n)]
        elif column == "guide_id":
            data[column] = [f"sg{i:04d}" for i in range(n)]
        elif column.endswith("_p_value"):
            data[column] = pvalues
        elif column.endswith("_negLog10p"):
            data[column] = np.full(n, np.nan)
        elif column.endswith(("_fc", "_se", "_value")):
            data[column] = rng.normal(size=n)
        elif column == "targeting":
            data[column] = rng.random(n) > 0.5
        elif column in {
            "guide_start", "guide_end", "element_start", "element_end",
            "intended_target_start", "intended_target_end",
            "nPerturbedCells", "num_guides",
        }:
            data[column] = rng.integers(1, 10**6, n)
        else:
            data[column] = [f"{column}_{i}" for i in range(n)]
    frame = pd.DataFrame(data)
    # mergedResults derives negLog10p with this helper before writing the table.
    for prefix in ("sceptre", "perturbo"):
        if f"{prefix}_p_value" in frame and f"{prefix}_negLog10p" in frame:
            frame[f"{prefix}_negLog10p"] = neg_log10(frame[f"{prefix}_p_value"], PVALUE_FLOOR)
    return frame


def _inputs(tmp_path, kwargs, *, with_local_perturbo=True):
    join = list(kwargs["join_columns"])
    n = 24
    pvalues = np.concatenate([np.linspace(1e-8, 0.9, n - 2), [1e-320, 0.0]])

    local_columns = join + list(kwargs["local_metric_columns"])
    if with_local_perturbo:
        local_columns += [c for c in LOCAL_PERTURBO_COLUMNS if c not in local_columns]
    local = _frame(local_columns, n, pvalues)
    if with_local_perturbo:
        local["perturbo_q_value"] = np.clip(local["perturbo_p_value"] * 1.5, 0.0, 1.0)
    glob = _frame(join + list(kwargs["global_required_columns"]), n, pvalues)

    local_path = tmp_path / "local.parquet"
    global_path = tmp_path / "global.parquet"
    local.to_parquet(local_path)
    glob.to_parquet(global_path)
    return local, local_path, global_path, tmp_path / "out.parquet"


@pytest.mark.parametrize("module,builder_attr", BUILDERS)
def test_fast_path_is_taken_when_local_has_perturbo_columns(module, builder_attr, tmp_path):
    kwargs = _fast_path_kwargs(module, builder_attr)
    local, local_path, global_path, out_path = _inputs(tmp_path, kwargs)

    taken = streaming_catalog_io.try_write_enriched_parquet_catalog(
        local_path=local_path,
        global_path=global_path,
        output_path=out_path,
        **{k: v for k, v in kwargs.items()
           if k not in {"local_path", "global_path", "output_path"}},
    )
    assert taken is True, "the catalog fell back to pandas"

    written = pd.read_parquet(out_path)
    assert list(written.columns) == list(kwargs["output_columns"])

    # The cis block must equal the local table's perturbo_* columns exactly --
    # it is a rename, so anything else means a value was recomputed or lost.
    keys = [c for c in kwargs["join_columns"] if c in written.columns]
    merged = written.merge(local, on=keys, suffixes=("", "_local"))
    for source, target in zip(LOCAL_PERTURBO_COLUMNS,
                              [f"perturbo_cis_{c[len('perturbo_'):]}" for c in LOCAL_PERTURBO_COLUMNS]):
        assert target in written.columns, target
        assert written[target].notna().any(), f"{target} came through empty"
        got, expected = merged[target], merged[f"{source}_local"]
        assert ((got == expected) | (got.isna() & expected.isna())).all(), target


@pytest.mark.parametrize("module,builder_attr", BUILDERS)
def test_fast_path_nulls_cis_block_without_local_perturbo(module, builder_attr, tmp_path):
    """A SCEPTRE-only run still takes the fast path; pandas fills NaN there."""
    kwargs = _fast_path_kwargs(module, builder_attr)
    _, local_path, global_path, out_path = _inputs(tmp_path, kwargs, with_local_perturbo=False)

    taken = streaming_catalog_io.try_write_enriched_parquet_catalog(
        local_path=local_path,
        global_path=global_path,
        output_path=out_path,
        **{k: v for k, v in kwargs.items()
           if k not in {"local_path", "global_path", "output_path"}},
    )
    assert taken is True
    written = pd.read_parquet(out_path)
    assert list(written.columns) == list(kwargs["output_columns"])
    cis_columns = [c for c in kwargs["output_columns"] if c.startswith("perturbo_cis_")]
    assert cis_columns
    for target in cis_columns:
        assert written[target].isna().all(), target


@pytest.mark.parametrize("module,builder_attr", BUILDERS)
def test_fast_path_defers_when_a_q_value_still_needs_filling(module, builder_attr, tmp_path):
    """Benjamini-Hochberg does not belong in a streaming plan, so defer instead."""
    kwargs = _fast_path_kwargs(module, builder_attr)
    local, local_path, global_path, out_path = _inputs(tmp_path, kwargs)
    local.loc[2, "perturbo_q_value"] = np.nan
    local.to_parquet(local_path)

    taken = streaming_catalog_io.try_write_enriched_parquet_catalog(
        local_path=local_path,
        global_path=global_path,
        output_path=out_path,
        **{k: v for k, v in kwargs.items()
           if k not in {"local_path", "global_path", "output_path"}},
    )
    assert taken is False


def test_collect_engine_keyword_matches_the_installed_polars():
    """A collect() keyword this Polars rejects would demote the whole catalog."""
    import inspect

    import polars as pl
    from polars_compat import COLLECT_ENGINE

    parameters = inspect.signature(pl.LazyFrame.collect).parameters
    for name in COLLECT_ENGINE:
        assert name in parameters, f"collect() has no {name!r} on polars {pl.__version__}"
