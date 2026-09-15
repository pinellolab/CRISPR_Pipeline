#!/usr/bin/env python
"""Score direct-target rows against matched non-targeting controls.

The evaluation reads one result table out of the inference MuData's ``uns``.
On the TAP-seq chr8 screen that table is 52,006,760 rows by 21 columns -- 4.99
GB of a 6.41 GB file -- and a handful of columns decide which rows the curves
are built from: ``guide_id`` carries the guide's intended target, ``gene_id``
and ``gene_name`` say whether the row is that target, and two metrics score it.
Fewer than a thousand rows survive.

So nothing here loads the table. ``mudata.read_h5mu`` is not used at all --
``backed`` defers only ``.X``, so it would read every result table in ``uns``
whatever this script goes on to touch. The three matching keys are read at
category resolution instead: the identifiers are stored as categories plus
codes, so the guide-to-target map and the identifier matching run over the
thousands of distinct names rather than over 52 million rows, and the full
columns are read only for the rows that reach a plot. The counts, the metrics
and the plots are unchanged.
"""

import mudata as md
from sklearn.metrics import precision_recall_curve, roc_curve, auc
import h5py
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

from inference_target_matching import direct_target_mask
from qc_mudata_io import read_mudata_without_uns, result_keys, result_table_columns


RESULTS_KEY = "global_analysis_per_guide_results"

# Above this many selected rows, reading a column in order and masking it beats
# asking HDF5 for that many scattered points.
_POINT_SELECT_LIMIT = 200_000
_SLAB_ROWS = 4_000_000


def select_inference_columns(df):
    """Return log2 fold-change and p-value columns available in inference output.

    Takes either a frame or its column names, because the names are all this
    needs and reading 52 million rows to look at them is the thing to avoid.
    """
    columns = list(df.columns) if hasattr(df, "columns") else list(df)
    candidate_pairs = [
        ("perturbo_log2_fc", "perturbo_p_value"),
        ("sceptre_log2_fc", "sceptre_p_value"),
        ("log2_fc", "p_value"),
    ]
    for fc_col, p_col in candidate_pairs:
        if fc_col in columns and p_col in columns:
            print(f"Using inference columns: {fc_col}, {p_col}")
            return fc_col, p_col

    raise KeyError(
        "Could not find inference columns. Expected one of: "
        + ", ".join(f"{fc}/{p}" for fc, p in candidate_pairs)
        + f". Available columns: {columns}"
    )


def savefig(path):
    """Utility to save figures cleanly."""
    plt.savefig(path, dpi=300, bbox_inches="tight")
    print(f"✔ Saved plot:", path)


def keep_scorable_rows(df, p_col):
    """Rows the curves can score, and how many were dropped.

    A pair the conditional randomization test did not test carries no p-value.
    Such a row cannot be scored, and dropping it after the classes have been
    matched below would move the precision-recall baseline off the intended
    1:1 prevalence, so both classes are filtered before they are counted.
    """
    values = pd.to_numeric(df[p_col], errors="coerce").astype(float)
    finite = np.isfinite(values.to_numpy())
    return df.loc[finite].copy(), int((~finite).sum())


def perform_binary_evaluation(label_controls, infered_significance_col, outdir, plot=True, evaluation_tag=''):
    evaluation_input = pd.DataFrame(
        {
            "label": label_controls,
            "p_value": pd.to_numeric(infered_significance_col, errors="coerce"),
        }
    )
    scorable = evaluation_input.dropna()
    unscorable = len(evaluation_input) - len(scorable)
    if unscorable:
        print(
            f"Dropped {unscorable} row(s) with a missing label or p-value before scoring."
        )
    evaluation_input = scorable
    true_label = evaluation_input["label"].astype(int)
    pred_value = evaluation_input["p_value"]

    if evaluation_input.empty or true_label.nunique() < 2:
        reason = (
            "Controls evaluation skipped: precision-recall and ROC curves require "
            "at least one finite p-value for both direct-target and control rows."
        )
        print(reason)
        with open(
            os.path.join(outdir, "controls_evaluation_skipped.txt"),
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(reason + "\n")
            handle.write(f"valid_rows={len(evaluation_input)}\n")
            handle.write(f"classes={sorted(true_label.unique().tolist())}\n")
        return False

    # 1 - p is the score, so a small p-value ranks as more likely a true
    # positive. The column must therefore be a p-value, never a q-value.
    pre, rec, _ = precision_recall_curve(true_label, 1-pred_value)
    auprc = auc(rec, pre)

    fpr, tpr, _ = roc_curve(true_label, 1-pred_value)
    auroc = auc(fpr, tpr)

    n_positive = int(true_label.sum())
    n_negative = int(len(true_label) - n_positive)
    positive_rate = n_positive / len(true_label)
    print(
        f"Scored {len(true_label)} rows: {n_positive} direct-target, "
        f"{n_negative} non-targeting (positive rate {positive_rate:.3f})"
    )
    print(f"Area under Precision-Recall Curve : {auprc:.3f}")
    print(f"Area under Receiver-Operating Curve: {auroc:.3f}")

    # The precision-recall baseline is the positive rate, so the class counts
    # the curves were actually built from belong beside the areas.
    with open(
        os.path.join(outdir, "controls_evaluation_summary.txt"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"direct_target_rows={n_positive}\n")
        handle.write(f"non_targeting_rows={n_negative}\n")
        handle.write(f"unscorable_rows_dropped={unscorable}\n")
        handle.write(f"positive_rate={positive_rate:.6f}\n")
        handle.write(f"auprc={auprc:.6f}\n")
        handle.write(f"auroc={auroc:.6f}\n")

    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=130)
        axes[0].plot(rec, pre, lw=1, label=f"AUPRC={auprc:.3f}")
        axes[1].plot(fpr, tpr, lw=1, label=f"AUROC={auroc:.3f}")

        axes[0].set_xlabel("Recall")
        axes[0].set_ylabel("Precision")
        axes[1].set_xlabel("False Positive Rate")
        axes[1].set_ylabel("True Positive Rate")
        axes[0].legend()
        axes[1].legend()
        plt.title(evaluation_tag)
        plt.tight_layout()

        # Save
        savefig(os.path.join(outdir, "global_analysis_perturbo_precision_recall_roc.png"))
        plt.show()
    return True


def plot_volcano(
    table_to_fdr,
    outdir,
    fc_col="log2_fc",
    p_col="p_value",
    direct_col="direct_target",
    figsize=(6, 8),
    fc_threshold=1,
    p_threshold=0.05,
    dpi=180
):
    df = table_to_fdr.copy()
    df["neglog10p"] = -np.log10(df[p_col].replace(0, np.nan))

    df["color"] = df[direct_col].map({1: "#d62728", 0: "#7f7f7f"})

    plt.figure(figsize=figsize, dpi=dpi)
    plt.scatter(
        df[fc_col],
        df["neglog10p"],
        c=df["color"],
        s=28,
        alpha=0.75,
        edgecolors="none"
    )

    plt.axvline(x=fc_threshold, color="black", linestyle="--", lw=1, alpha=0.6)
    plt.axvline(x=-fc_threshold, color="black", linestyle="--", lw=1, alpha=0.6)
    plt.axhline(y=-np.log10(p_threshold), color="black", linestyle="--", lw=1, alpha=0.6)

    plt.xlabel("log2 Fold Change", fontsize=13)
    plt.ylabel("-log10(p-value)", fontsize=13)
    plt.title("Volcano Plot", fontsize=15, weight="bold")

    plt.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.35)
    plt.gca().set_facecolor("#f7f7f7")
    for spine in ["top", "right"]:
        plt.gca().spines[spine].set_visible(False)

    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color="#d62728", label="Guides Direct target")
    gray_patch = mpatches.Patch(color="#7f7f7f", label="Non-targeting")
    plt.legend(handles=[red_patch, gray_patch], frameon=False)

    plt.tight_layout()

    # Save
    savefig(os.path.join(outdir, "global_analysis_perturbo_volcano_plot.png"))
    plt.show()


def _decode(values):
    """Decode the byte strings h5py hands back for a text dataset."""
    if values.dtype.kind in "SO":
        return np.array(
            [v.decode() if isinstance(v, bytes) else v for v in values], dtype=object
        )
    return values


def _encoding(node):
    encoding = node.attrs.get("encoding-type", b"")
    return encoding.decode() if isinstance(encoding, bytes) else str(encoding)


def _codes_dtype(n_categories):
    """The narrowest integer type that can hold these codes and -1."""
    for dtype in (np.int8, np.int16, np.int32):
        if n_categories <= np.iinfo(dtype).max:
            return dtype
    return np.int64


def _has_missing(codes):
    return bool(codes.size) and bool(codes.min() < 0)


def _is_missing(value):
    return value is None or (isinstance(value, float) and np.isnan(value))


def _equals_false(value):
    """``value == False``, the comparison the whole-column query would make."""
    return bool(pd.Series([value], dtype=object).eq(False).iloc[0])


def _row_values(per_category, codes, na_value):
    """One value per row, from one value per category.

    The whole point of keeping the identifiers categorical: the work above
    happens once per distinct name, and only this gather is per row.
    """
    if not len(per_category):
        return np.full(len(codes), na_value, dtype=np.asarray([na_value]).dtype)
    values = per_category[np.where(codes >= 0, codes, 0)]
    missing = codes < 0
    if missing.any():
        values[missing] = na_value
    return values


def _mapped_categories(categories, mapping, has_missing_code):
    """``Series.map`` applied to the categories, the way pandas applies it.

    ``Series.map`` on a categorical column maps the categories too, but it
    materializes one object per row whenever the mapping is not one-to-one,
    which for 52 million rows is the expensive part of an otherwise tiny
    lookup. ``Index.map`` performs the same lookup -- missing keys become
    missing values -- and the missing category takes the value pandas gives it.
    """
    mapped = np.asarray(pd.Index(categories).map(mapping), dtype=object)
    na_value = mapping.get(np.nan, np.nan) if has_missing_code else np.nan
    return mapped, na_value


def _map_categorical(values, mapping):
    """``Series.map(mapping)`` over a categorical column, without leaving it."""
    codes = values.codes
    mapped, na_value = _mapped_categories(
        values.categories, mapping, _has_missing(codes)
    )
    inner, uniques = pd.factorize(mapped, sort=False)
    uniques = np.asarray(uniques, dtype=object)
    if _is_missing(na_value):
        na_code = -1
    else:
        found = np.flatnonzero(uniques == na_value)
        if len(found):
            na_code = int(found[0])
        else:
            uniques = np.append(uniques, na_value)
            na_code = len(uniques) - 1
    inner = inner.astype(_codes_dtype(len(uniques)), copy=False)
    return pd.Categorical.from_codes(
        _row_values(inner, codes, na_code), categories=pd.Index(uniques)
    )


def _as_categorical(values):
    """A column as a Categorical, whatever dtype it arrived in.

    ``factorize(sort=False)`` never compares two values, so a column holding
    mixed types still resolves the way the string comparison below would.
    """
    if isinstance(values.dtype, pd.CategoricalDtype):
        return pd.Categorical(values)
    codes, uniques = pd.factorize(np.asarray(values, dtype=object), sort=False)
    return pd.Categorical.from_codes(
        codes.astype(_codes_dtype(len(uniques)), copy=False),
        categories=pd.Index(np.asarray(uniques, dtype=object)),
    )


class _FrameResults:
    """A result table that is already in memory."""

    def __init__(self, frame):
        self._frame = frame

    @property
    def columns(self):
        return list(self._frame.columns)

    def __len__(self):
        return len(self._frame)

    def key_column(self, name):
        return _as_categorical(self._frame[name])

    def take(self, positions, columns):
        return self._frame.iloc[positions][list(columns)]


class _H5Results:
    """A result table left on disk, read one column and one row set at a time.

    ``qc_mudata_io`` supplies the names and the modalities; the reads below are
    here because they are by row set, which the shared reader does not do -- it
    reads whole columns, and a whole column is what this is avoiding.
    """

    def __init__(self, path, key):
        self._path = path
        self._key = key
        self._columns = result_table_columns(path, key)
        with h5py.File(path, "r") as handle:
            group = handle[f"uns/{key}"]
            self._n_rows = self._length(group)

    @staticmethod
    def _length(group):
        for name in ("_index", *group.keys()):
            node = group.get(name)
            if isinstance(node, h5py.Group):
                node = node.get("codes", node.get("values"))
            if node is not None:
                return int(node.shape[0])
        return 0

    @property
    def columns(self):
        return list(self._columns)

    def __len__(self):
        return self._n_rows

    def key_column(self, name):
        with h5py.File(self._path, "r") as handle:
            node = handle[f"uns/{self._key}"][name]
            if _encoding(node) == "categorical":
                # Categories plus codes is exactly the representation the
                # matching wants, so it is read as it is stored.
                categories = pd.Index(_decode(node["categories"][:]))
                return pd.Categorical.from_codes(node["codes"][:], categories=categories)
            return _as_categorical(pd.Series(self._whole(node)))

    def take(self, positions, columns):
        positions = np.asarray(positions, dtype=np.int64)
        mask = None
        if len(positions) > _POINT_SELECT_LIMIT:
            mask = np.zeros(self._n_rows, dtype=bool)
            mask[positions] = True
        with h5py.File(self._path, "r") as handle:
            group = handle[f"uns/{self._key}"]
            frame = {
                name: self._rows(group[name], positions, mask) for name in columns
            }
            index = group.get("_index")
            index = (
                _decode(self._read(index, positions, mask))
                if index is not None
                else positions
            )
        return pd.DataFrame(frame, index=pd.Index(index))

    def _read(self, node, positions, mask):
        """One column's values at the selected rows."""
        if mask is None:
            return node[positions] if len(positions) else node[0:0]
        blocks = [
            node[start : start + _SLAB_ROWS][mask[start : start + _SLAB_ROWS]]
            for start in range(0, self._n_rows, _SLAB_ROWS)
        ]
        return np.concatenate(blocks) if blocks else node[0:0]

    def _rows(self, node, positions, mask):
        """One column at the selected rows, decoded as anndata wrote it."""
        encoding = _encoding(node)
        if encoding == "categorical":
            categories = _decode(node["categories"][:])
            codes = self._read(node["codes"], positions, mask)
            values = np.full(len(codes), None, dtype=object)
            present = codes >= 0
            values[present] = categories[codes[present]]
            return values
        if encoding.startswith("nullable"):
            return self._nullable(
                self._read(node["values"], positions, mask),
                self._read(node["mask"], positions, mask),
            )
        return _decode(self._read(node, positions, mask))

    def _whole(self, node):
        return self._rows(node, None, np.ones(self._n_rows, dtype=bool))

    @staticmethod
    def _nullable(values, mask):
        """Values plus mask, back to the nullable dtype they were written as."""
        mask = mask.astype(bool, copy=False)
        if values.dtype.kind in "iu":
            return pd.arrays.IntegerArray(values, mask)
        if values.dtype.kind == "b":
            return pd.arrays.BooleanArray(values, mask)
        if values.dtype.kind == "f":
            values[mask] = np.nan
            return values
        values = values.astype(object)
        values[mask] = None
        return values


def _write_skip(outdir, reason, *lines):
    print(reason)
    with open(
        os.path.join(outdir, "controls_evaluation_skipped.txt"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(reason + "\n")
        for line in lines:
            handle.write(line + "\n")


def run_evaluation_controls(md_read, outdir):
    """Evaluate the controls from an already-loaded MuData."""
    os.makedirs(outdir, exist_ok=True)

    col_used = RESULTS_KEY
    if col_used not in md_read.uns:
        _write_skip(
            outdir,
            "Controls evaluation skipped: global PerTurbo guide results are "
            "not present for this inference mode.",
            "global_analysis_per_guide_results=absent",
        )
        return
    evaluate_controls_table(
        _FrameResults(md_read.uns[col_used]), md_read['guide'].var, outdir
    )


def run_evaluation_controls_from_path(mdata_path, outdir):
    """Evaluate the controls without loading the result tables out of ``uns``."""
    os.makedirs(outdir, exist_ok=True)

    print("Loading MuData file...")
    md_read = read_mudata_without_uns(mdata_path)
    available = result_keys(mdata_path)
    print("Finished Loading MuData file...")

    if RESULTS_KEY not in available:
        _write_skip(
            outdir,
            "Controls evaluation skipped: global PerTurbo guide results are "
            "not present for this inference mode.",
            "global_analysis_per_guide_results=absent",
        )
        return
    evaluate_controls_table(
        _H5Results(mdata_path, RESULTS_KEY), md_read['guide'].var, outdir
    )


def evaluate_controls_table(results, guide_var, outdir):
    """Score the direct-target rows of one result table against its controls."""
    fc_col, p_col = select_inference_columns(results.columns)
    #converting to avoid non boolean values
    col = guide_var['targeting']

    guide_var['targeting'] = col.apply(
        lambda x: True if (x is True or str(x).upper() == "TRUE")
        else False if (x is False or str(x).upper() == "FALSE")
        else x
    )

    selecting_non_targeting_guides = guide_var[guide_var['targeting'] == False]
    non_targeting_ids = set(selecting_non_targeting_guides['guide_id'].values)
    print(f"Number of non-targeting guides: {len(non_targeting_ids)}")
    intended_dict = guide_var.set_index(['guide_id'])['intended_target_name'].to_dict()
    print(f"Number of unique intended targets: {len(set(intended_dict.values()))}")

    columns = results.columns
    guide_ids = results.key_column('guide_id')
    # The guide's intended target replaces whatever the table carries, as it
    # always has; the lookup just happens once per guide instead of once per row.
    intended_target_name = _map_categorical(guide_ids, intended_dict)

    keys = {'intended_target_name': intended_target_name}
    for name in ('gene_id', 'gene_name'):
        if name in columns:
            keys[name] = results.key_column(name)
    direct_mask = direct_target_mask(pd.DataFrame(keys)).to_numpy()
    direct_positions = np.flatnonzero(direct_mask)

    selecting_to_plot_final = _selected_rows(
        results, direct_positions, columns, intended_target_name
    ).drop_duplicates()

    dict_targeting_or_no = guide_var.set_index('guide_id')['targeting'].to_dict()
    targeting_genes, targeting_na = _mapped_categories(
        guide_ids.categories, dict_targeting_or_no, _has_missing(guide_ids.codes)
    )
    # "== False" per distinct guide, then per row: the same comparison the
    # whole-column query made, over 4,120 values instead of 52 million.
    is_control = (pd.Series(targeting_genes, dtype=object) == False).to_numpy()
    control_mask = _row_values(
        is_control, guide_ids.codes, _equals_false(targeting_na)
    )

    gene_ids = keys['gene_id']
    all_targets = selecting_to_plot_final['gene_id'].dropna().drop_duplicates()
    print (f"all targets: {all_targets.shape[0]}")
    using_random_sample = ''
    matching_mask = control_mask & _row_values(
        pd.Index(gene_ids.categories).isin(all_targets), gene_ids.codes, False
    )
    if matching_mask.any():
        control_positions = np.flatnonzero(matching_mask)
    else:
        control_positions = np.flatnonzero(control_mask)
        if len(control_positions) > 0:
            print ("No non-targeting control guides found for evaluation that were tested against the intended targets.")
            print ('Using random sample of non-targeting guides for evaluation instead.')
            using_random_sample = '\n Warning \n Using random sample of non-targeting guides for evaluation instead.'

    # A control row reaches the output only through its guide, its fold change
    # and its p-value, so those are the only columns read for it. A plot that
    # starts reading some other column has to ask for it here.
    control_columns = [
        name for name in ('guide_id', 'gene_id', fc_col, p_col) if name in columns
    ]
    non_target_controls = results.take(control_positions, control_columns).copy()
    non_target_controls['targeting_genes'] = _row_values(
        targeting_genes, guide_ids.codes[control_positions], targeting_na
    )
    non_target_controls['direct_target'] = 0
    print (f"Number of non-targeting control guides selected for evaluation: {non_target_controls.shape[0]}")

    table_to_test_cis = selecting_to_plot_final.copy()
    print (_selected_rows(
        results, np.arange(min(5, len(results))), columns, intended_target_name
    ).values)

    table_to_test_cis['direct_target'] = 1
    print (f"Number of targeting guides for direct targets: {table_to_test_cis.shape[0]}")

    # Both classes are cut down to the rows that carry a finite p-value before
    # they are counted or matched. Since the conditional randomization test
    # stopped passing the posterior probability off as a p-value, an untested
    # pair is missing rather than filled, and matching the classes on counts
    # that include unscorable rows would leave the curves unbalanced.
    table_to_test_cis, dropped_targets = keep_scorable_rows(table_to_test_cis, p_col)
    non_target_controls, dropped_controls = keep_scorable_rows(non_target_controls, p_col)
    if dropped_targets or dropped_controls:
        print(
            f"Rows without a finite {p_col} dropped before matching: "
            f"{dropped_targets} direct-target, {dropped_controls} non-targeting."
        )
        print(
            f"Scorable rows: {table_to_test_cis.shape[0]} direct-target, "
            f"{non_target_controls.shape[0]} non-targeting."
        )

    if table_to_test_cis.empty or non_target_controls.empty:
        missing_group = "direct-target rows" if table_to_test_cis.empty else "non-targeting guides"
        _write_skip(
            outdir,
            f"Controls evaluation skipped: no {missing_group} with a finite "
            f"{p_col} are present in the inference results, so AUROC/AUPRC and "
            "matched-control plots cannot be calculated.",
            f"targeting_direct_target_rows={table_to_test_cis.shape[0]}",
            f"non_targeting_control_rows={non_target_controls.shape[0]}",
        )
        return

    sample_with_replacement = non_target_controls.shape[0] < table_to_test_cis.shape[0]
    if sample_with_replacement:
        print(
            "Fewer non-targeting controls than direct-target rows; "
            "sampling controls with replacement."
        )
    table_to_fdr = pd.concat([
        table_to_test_cis,
        non_target_controls.sample(
            n=table_to_test_cis.shape[0],
            random_state=42,
            replace=sample_with_replacement,
        )
    ]).copy()
    table_to_fdr["log2_fc"] = table_to_fdr[fc_col]
    table_to_fdr["p_value"] = table_to_fdr[p_col]
    print (table_to_fdr)
    # Volcano plot
    plot_volcano(table_to_fdr, outdir=outdir)

    # Binary evaluation curves
    perform_binary_evaluation(
        table_to_fdr['direct_target'],
        table_to_fdr['p_value'],
        outdir=outdir,
        plot=True,
        evaluation_tag = using_random_sample,
    )

    # Bar plot
    plt.figure(figsize=(5, 4), dpi=150)
    table_to_fdr.groupby('direct_target').count()['guide_id'].plot(kind='bar')
    plt.ylabel('Number of guides')
    plt.title(f"Direct targets vs random control guides \n {using_random_sample}")
    plt.tight_layout()

    savefig(os.path.join(outdir, "global_analysis_perturbo_barplot_direct_vs_control.png"))
    plt.show()


def _selected_rows(results, positions, columns, intended_target_name):
    """The full row for each selected position, intended target filled in."""
    frame = results.take(positions, columns).copy()
    frame['intended_target_name'] = np.asarray(intended_target_name.take(positions))
    return frame


if __name__ == "__main__":
    print("running controls evaluation program...")
    parser = argparse.ArgumentParser(description="Running controls evaluation program")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--outdir", type=str, default="plots", help="Directory to save plots")

    args = parser.parse_args()

    # The result tables in uns are 6.35 GB of the 6.41 GB TAP-seq chr8 file and
    # read_h5mu loads all of them -- backed defers only .X -- so the modalities
    # are read on their own and the one result table is read column by column.
    run_evaluation_controls_from_path(args.mdata_path, outdir=args.outdir)
