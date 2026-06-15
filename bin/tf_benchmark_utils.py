import pandas as pd


def normalize_tf_benchmark_results(results: pd.DataFrame) -> pd.DataFrame:
    """Normalize supported trans-result schemas for the TF benchmark."""
    normalized = results.copy()

    aliases = {
        "p_value": ["perturbo_p_value"],
        "log2_fc": ["perturbo_log2_fc"],
    }
    for canonical, candidates in aliases.items():
        for candidate in candidates:
            if candidate not in normalized.columns:
                continue
            if canonical not in normalized.columns:
                normalized[canonical] = normalized[candidate]
            else:
                normalized[canonical] = normalized[canonical].fillna(
                    normalized[candidate]
                )

    if "p_value" not in normalized.columns:
        raise ValueError(
            "TF benchmark requires a p-value column. Expected one of: "
            "p_value, perturbo_p_value. "
            f"Available columns: {sorted(normalized.columns.tolist())}"
        )

    normalized["p_value"] = pd.to_numeric(
        normalized["p_value"], errors="coerce"
    )
    if not normalized["p_value"].notna().any():
        raise ValueError("TF benchmark p-value column contains no numeric values.")

    return normalized
