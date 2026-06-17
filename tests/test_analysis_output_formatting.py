import sys
from pathlib import Path

import numpy as np
import pandas as pd


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from analysis_output_formatting import add_neg_log10_columns, make_h5mu_safe_dataframe


def test_add_neg_log10_columns_for_method_specific_pvalues():
    results = pd.DataFrame(
        {
            "sceptre_p_value": [0.01, np.nan],
            "perturbo_p_value": [0.2, 0.0],
        }
    )

    formatted = add_neg_log10_columns(results)

    assert np.isclose(formatted.loc[0, "sceptre_negLog10p"], 2.0)
    assert formatted.loc[0, "sceptre_log10_p_value"] == formatted.loc[0, "sceptre_negLog10p"]
    assert pd.isna(formatted.loc[1, "sceptre_negLog10p"])
    assert np.isclose(formatted.loc[0, "perturbo_negLog10p"], -np.log10(0.2))
    assert formatted.loc[1, "perturbo_negLog10p"] == 300.0
    assert formatted.loc[1, "perturbo_log10_p_value"] == 300.0


def test_make_h5mu_safe_dataframe_converts_nullable_strings_to_object():
    results = pd.DataFrame(
        {
            "intended_target_name": pd.Series(["target1", pd.NA], dtype="string"),
            "p_value": [0.01, 0.2],
        }
    )

    safe = make_h5mu_safe_dataframe(results)

    assert safe["intended_target_name"].dtype == object
    assert safe.loc[0, "intended_target_name"] == "target1"
    assert safe.loc[1, "intended_target_name"] is None
