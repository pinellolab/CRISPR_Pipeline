import numpy as np
import pandas as pd
from scipy import sparse


ELEMENT_COLUMNS = [
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
]
ELEMENT_ANNOTATION_COLUMNS = [
    "element_id",
    "element_type",
    "element_chr",
    "element_start",
    "element_end",
    "element_name",
    "guide_ids",
    "num_guides",
    "gene_name",
    "nPerturbedCells",
]
SCEPTRE_RESULT_COLUMNS = [
    "sceptre_log2_fc",
    "sceptre_p_value",
    "sceptre_q_value",
    "sceptre_fc_se",
    "sceptre_negLog10p",
]
PERTURBO_RESULT_COLUMNS = [
    "perturbo_log2_fc",
    "perturbo_p_value",
    "perturbo_q_value",
    "perturbo_fc_se",
    "perturbo_negLog10p",
]
GUIDE_ANNOTATION_COLUMNS = [
    "guide_sequence",
    "guide_type",
    "targeting",
    "guide_chr",
    "guide_start",
    "guide_end",
    "guide_strand",
    "pam",
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
    "gene_name",
    "nPerturbedCells",
]
P_VALUE_FLOOR = 1e-300
MISSING_TOKEN = ""


def _normalize_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def _normalize_int(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def _normalize_element_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ELEMENT_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA

    out["intended_target_name"] = _normalize_text(out["intended_target_name"])
    out["intended_target_chr"] = _normalize_text(out["intended_target_chr"])
    out["intended_target_start"] = _normalize_int(out["intended_target_start"])
    out["intended_target_end"] = _normalize_int(out["intended_target_end"])
    return out


def _element_key(df: pd.DataFrame) -> pd.Series:
    normalized = _normalize_element_columns(df)
    key_parts = [
        normalized["intended_target_name"].fillna(MISSING_TOKEN),
        normalized["intended_target_chr"].fillna(MISSING_TOKEN),
        normalized["intended_target_start"].astype("string").fillna(MISSING_TOKEN),
        normalized["intended_target_end"].astype("string").fillna(MISSING_TOKEN),
    ]
    return key_parts[0] + "|" + key_parts[1] + "|" + key_parts[2] + "|" + key_parts[3]


def _collapse_unique_strings(values: pd.Series) -> object:
    cleaned = sorted(
        {str(v).strip() for v in values if pd.notna(v) and str(v).strip() != ""}
    )
    if not cleaned:
        return pd.NA
    return ";".join(cleaned)


def neg_log10(series: pd.Series, pvalue_floor: float = P_VALUE_FLOOR) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    clamped = values.where(values.notna()).clip(lower=pvalue_floor)
    return -np.log10(clamped)


def add_neg_log10_columns(
    df: pd.DataFrame, pvalue_floor: float = P_VALUE_FLOOR
) -> pd.DataFrame:
    out = df.copy()
    for prefix in ("sceptre", "perturbo"):
        pvalue_col = f"{prefix}_p_value"
        if pvalue_col not in out.columns:
            continue
        neglog_col = f"{prefix}_negLog10p"
        out[neglog_col] = neg_log10(out[pvalue_col], pvalue_floor)

    if (
        "p_value" in out.columns
        and "sceptre_p_value" not in out.columns
        and "perturbo_p_value" not in out.columns
    ):
        out["negLog10p"] = neg_log10(out["p_value"], pvalue_floor)

    return out


def make_h5mu_safe_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert pandas extension string columns before storing tables in MuData.uns."""
    out = df.copy()
    for col in out.columns:
        series = out[col]
        if pd.api.types.is_string_dtype(series.dtype) or pd.api.types.is_object_dtype(
            series.dtype
        ):
            out[col] = series.astype(object).where(series.notna(), None)
    return out


def _get_modality(mdata, name: str):
    if hasattr(mdata, "mod") and name in mdata.mod:
        return mdata[name]
    return None


def _empty_element_metadata() -> pd.DataFrame:
    columns = [
        "_element_key",
        "guide_ids",
        "num_guides",
        "element_type",
        "nPerturbedCells",
    ]
    return pd.DataFrame(columns=columns)


def _empty_guide_metadata() -> pd.DataFrame:
    return pd.DataFrame(columns=["guide_id"] + GUIDE_ANNOTATION_COLUMNS[:-2] + ["nPerturbedCells"])


def _build_guide_metadata(mdata) -> pd.DataFrame:
    guide = _get_modality(mdata, "guide")
    if guide is None or guide.var.empty:
        return _empty_guide_metadata()

    guide_var = guide.var.copy()
    if "guide_id" not in guide_var.columns:
        guide_var["guide_id"] = guide_var.index.astype(str)
    guide_var["guide_id"] = guide_var["guide_id"].astype(str)
    if guide_var["guide_id"].duplicated().any():
        duplicates = guide_var.loc[
            guide_var["guide_id"].duplicated(keep=False), "guide_id"
        ].head(10).tolist()
        raise ValueError(f"guide metadata contains duplicate guide_id values: {duplicates}")

    aliases = {
        "guide_sequence": ["guide_sequence", "spacer", "sequence"],
        "guide_type": ["guide_type", "type"],
        "targeting": ["targeting"],
        "guide_chr": ["guide_chr"],
        "guide_start": ["guide_start"],
        "guide_end": ["guide_end"],
        "guide_strand": ["guide_strand", "strand"],
        "pam": ["pam"],
        "intended_target_name": ["intended_target_name"],
        "intended_target_chr": ["intended_target_chr"],
        "intended_target_start": ["intended_target_start"],
        "intended_target_end": ["intended_target_end"],
    }
    metadata = pd.DataFrame({"guide_id": guide_var["guide_id"].to_numpy()})
    numeric_columns = {
        "guide_start",
        "guide_end",
        "intended_target_start",
        "intended_target_end",
    }
    for output_col, candidates in aliases.items():
        source_col = next((col for col in candidates if col in guide_var.columns), None)
        if source_col is None:
            if output_col in numeric_columns:
                metadata[output_col] = np.nan
            elif output_col == "targeting":
                metadata[output_col] = False
            else:
                metadata[output_col] = ""
        elif output_col in numeric_columns:
            metadata[output_col] = pd.to_numeric(
                guide_var[source_col], errors="coerce"
            ).to_numpy()
        elif output_col == "targeting":
            metadata[output_col] = guide_var[source_col].map(
                lambda value: (
                    value
                    if isinstance(value, (bool, np.bool_))
                    else str(value).strip().lower() in {"true", "1", "t", "yes"}
                )
                if pd.notna(value)
                else False
            ).to_numpy(dtype=bool)
        else:
            metadata[output_col] = (
                guide_var[source_col].astype("string").fillna("").to_numpy(dtype=object)
            )

    assignment = guide.layers["guide_assignment"] if "guide_assignment" in guide.layers else guide.X
    if sparse.issparse(assignment):
        counts = np.asarray((assignment > 0).sum(axis=0)).ravel()
    else:
        counts = np.asarray(assignment > 0).sum(axis=0).ravel()
    if len(counts) != len(metadata):
        raise ValueError(
            "Guide assignment columns do not match the number of guide metadata rows."
        )
    metadata["nPerturbedCells"] = pd.array(counts, dtype="Int64")
    return metadata


def _build_element_metadata(mdata) -> pd.DataFrame:
    guide = _get_modality(mdata, "guide")
    if guide is None:
        return _empty_element_metadata()

    guide_var = guide.var.copy()
    if guide_var.empty:
        return _empty_element_metadata()
    if "guide_id" not in guide_var.columns:
        guide_var["guide_id"] = guide_var.index.astype(str)

    guide_var = _normalize_element_columns(guide_var)
    if "type" not in guide_var.columns:
        guide_var["type"] = pd.NA

    guide_var["guide_id"] = guide_var["guide_id"].astype("string")
    guide_var["_element_key"] = _element_key(guide_var)

    grouped = (
        guide_var.groupby(ELEMENT_COLUMNS + ["_element_key"], dropna=False, as_index=False)
        .agg(
            guide_ids=("guide_id", _collapse_unique_strings),
            num_guides=("guide_id", "nunique"),
            element_type=("type", _collapse_unique_strings),
        )
    )

    assignment = guide.layers["guide_assignment"] if "guide_assignment" in guide.layers else guide.X
    if not sparse.issparse(assignment):
        assignment = sparse.csr_matrix(assignment)
    else:
        assignment = assignment.tocsr()

    element_codes, unique_keys = pd.factorize(guide_var["_element_key"], sort=False)
    indicator = sparse.csr_matrix(
        (
            np.ones(len(element_codes), dtype=np.int8),
            (np.arange(len(element_codes)), element_codes),
        ),
        shape=(len(element_codes), len(unique_keys)),
    )
    cell_by_element = assignment @ indicator
    n_perturbed = np.asarray((cell_by_element > 0).sum(axis=0)).ravel().astype(int)
    n_perturbed_map = dict(zip(unique_keys.tolist(), n_perturbed.tolist()))

    grouped["nPerturbedCells"] = grouped["_element_key"].map(n_perturbed_map).astype("Int64")
    grouped["num_guides"] = grouped["num_guides"].astype("Int64")
    return grouped[
        ["_element_key", "guide_ids", "num_guides", "element_type", "nPerturbedCells"]
    ]


def _build_gene_name_map(mdata) -> pd.Series:
    for modality_name in ("gene", "rna"):
        modality = _get_modality(mdata, modality_name)
        if modality is None:
            continue
        gene_var = modality.var.copy()
        if gene_var.empty:
            continue

        if "symbol" in gene_var.columns:
            name_col = "symbol"
        elif "gene_name" in gene_var.columns:
            name_col = "gene_name"
        else:
            continue

        mapping = pd.Series(
            gene_var[name_col].astype("string").values,
            index=gene_var.index.astype(str),
        )
        return mapping[~mapping.index.duplicated(keep="first")]

    return pd.Series(dtype="string")


def _fill_gene_names(results: pd.DataFrame, mdata) -> pd.Series:
    gene_name_map = _build_gene_name_map(mdata)
    if gene_name_map.empty or "gene_id" not in results.columns:
        return pd.Series(pd.NA, index=results.index, dtype="string")

    out = results["gene_id"].astype(str).map(gene_name_map)
    if out.notna().all():
        return out

    stripped_index = gene_name_map.index.to_series().str.split(".").str[0]
    stripped_map = (
        pd.DataFrame({"gene_id": stripped_index.values, "gene_name": gene_name_map.values})
        .dropna(subset=["gene_id"])
        .drop_duplicates(subset=["gene_id"], keep="first")
        .set_index("gene_id")["gene_name"]
    )
    missing_mask = out.isna()
    out.loc[missing_mask] = (
        results.loc[missing_mask, "gene_id"]
        .astype(str)
        .str.split(".")
        .str[0]
        .map(stripped_map)
    )
    return out


def add_element_annotations(df: pd.DataFrame, mdata) -> pd.DataFrame:
    out = _normalize_element_columns(df).drop(
        columns=ELEMENT_ANNOTATION_COLUMNS + ["_element_key"],
        errors="ignore",
    )
    out["_element_key"] = _element_key(out)

    element_meta = _build_element_metadata(mdata)
    out = out.merge(element_meta, on="_element_key", how="left")

    out["element_name"] = out["intended_target_name"]
    out["element_id"] = out["element_name"]
    out["element_chr"] = out["intended_target_chr"]
    out["element_start"] = out["intended_target_start"]
    out["element_end"] = out["intended_target_end"]
    out["gene_name"] = _fill_gene_names(out, mdata)

    out = out.drop(columns=["_element_key"], errors="ignore")
    existing_annotations = [col for col in ELEMENT_ANNOTATION_COLUMNS if col in out.columns]
    non_annotations = [col for col in out.columns if col not in existing_annotations]
    return out[non_annotations + existing_annotations]


def add_guide_annotations(df: pd.DataFrame, mdata) -> pd.DataFrame:
    out = df.drop(columns=GUIDE_ANNOTATION_COLUMNS, errors="ignore").copy()
    out["guide_id"] = out["guide_id"].astype(str)
    out = out.merge(_build_guide_metadata(mdata), on="guide_id", how="left")
    out["gene_name"] = _fill_gene_names(out, mdata).fillna("")
    return out


def format_guide_output(df: pd.DataFrame, mdata=None) -> pd.DataFrame:
    out = add_neg_log10_columns(df)
    if mdata is not None:
        out = add_guide_annotations(out, mdata)
    preferred = (
        ["gene_id", "guide_id"]
        + SCEPTRE_RESULT_COLUMNS
        + PERTURBO_RESULT_COLUMNS
        + GUIDE_ANNOTATION_COLUMNS
    )
    ordered = [col for col in preferred if col in out.columns]
    extras = [col for col in out.columns if col not in ordered]
    return out[ordered + extras]


def format_element_output(df: pd.DataFrame, mdata) -> pd.DataFrame:
    out = add_element_annotations(add_neg_log10_columns(df), mdata)
    preferred = (
        ["gene_id"]
        + ELEMENT_COLUMNS
        + SCEPTRE_RESULT_COLUMNS
        + PERTURBO_RESULT_COLUMNS
        + ELEMENT_ANNOTATION_COLUMNS
    )
    ordered = [col for col in preferred if col in out.columns]
    extras = [col for col in out.columns if col not in ordered]
    return out[ordered + extras]
