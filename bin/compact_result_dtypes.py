"""Optional storage rounding for bounded result batches, never for inference.

Linear p/q values and posterior probabilities deliberately stay untouched.
Callers must compute statistics before applying this helper and preserve the
unrounded raw result artifacts. Unknown columns are also left unchanged.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# Explicit public fields only: suffix heuristics could accidentally cast a new
# probability or test statistic that requires different numerical guarantees.
FLOAT32_STORAGE_COLUMNS = frozenset(
    {"log2_fc", "perturbo_fc_se"}
    | {
        f"{prefix}_{field}"
        for prefix in ("perturbo", "perturbo_cis", "sceptre")
        for field in ("log2_fc", "fc_se", "negLog10p", "negLog10q")
    }
)


def compact_result_floats(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a batch with allowlisted storage fields rounded to float32.

    Other columns share their buffers; this function never mutates the input.
    Explicitly reject overflow/underflow of a finite nonzero number rather than
    silently replacing it with infinity or zero. NaNs and existing infinities
    retain their meaning. Apply only to bounded batches, after BH correction.
    """
    out = frame.copy(deep=False)
    for column in sorted(FLOAT32_STORAGE_COLUMNS.intersection(frame.columns)):
        values = pd.to_numeric(frame[column], errors="raise").to_numpy(
            dtype=np.float64, na_value=np.nan
        )
        with np.errstate(over="ignore", under="ignore"):
            packed = values.astype(np.float32)
        finite = np.isfinite(values)
        if np.any(finite & ~np.isfinite(packed)) or np.any(
            finite & (values != 0) & (packed == 0)
        ):
            raise ValueError(
                f"Cannot store {column} as float32 without overflow or underflow; "
                "disable compact result floats."
            )
        out[column] = packed
    return out
