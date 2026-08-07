import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "bin" / "intended_target_key_utils.py"
SPEC = importlib.util.spec_from_file_location("intended_target_key_utils", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONCAT_MODULE_PATH = Path(__file__).parents[1] / "bin" / "mudata_concat.py"
CONCAT_SPEC = importlib.util.spec_from_file_location("mudata_concat", CONCAT_MODULE_PATH)
CONCAT_MODULE = importlib.util.module_from_spec(CONCAT_SPEC)
CONCAT_SPEC.loader.exec_module(CONCAT_MODULE)


def test_explicit_control_element_ids_preserve_dual_guide_pairs():
    guide_var = pd.DataFrame(
        {
            "guide_id": ["target_a", "target_b", "nt_z", "nt_a", "nt_y", "nt_b"],
            "targeting": [True, True, False, False, False, False],
            "type": ["targeting", "targeting", "non_targeting_control", "non_targeting_control", "non_targeting_control", "non_targeting_control"],
            "intended_target_name": ["GENE1", "GENE1", "non-targeting", "non-targeting", "non-targeting", "non-targeting"],
            "intended_target_chr": ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
            "intended_target_start": [1, 1, 1, 1, 1, 1],
            "intended_target_end": [2, 2, 2, 2, 2, 2],
            "element_id": ["target_element", "target_element", "control_two", "control_one", "control_two", "control_one"],
        }
    )

    result = MODULE.annotate_intended_target_groups(guide_var)
    control_names = result.set_index("guide_id")["intended_target_name"]

    assert control_names["nt_z"] == control_names["nt_y"]
    assert control_names["nt_a"] == control_names["nt_b"]
    assert control_names["nt_z"] != control_names["nt_a"]


def test_concat_restores_source_only_element_id_annotation():
    source = pd.DataFrame(
        {
            "guide_id": ["ctrl_a", "ctrl_b"],
            "element_id": ["source_pair_1", "source_pair_1"],
        },
        index=["ctrl_a", "ctrl_b"],
    )
    combined = pd.DataFrame(
        {"guide_id": ["ctrl_a", "ctrl_b"]}, index=["ctrl_a", "ctrl_b"]
    )

    restored = CONCAT_MODULE.preserve_source_guide_metadata(combined, source)

    assert restored["element_id"].tolist() == ["source_pair_1", "source_pair_1"]


def test_concat_does_not_coerce_existing_annotation_dtype():
    source = pd.DataFrame(
        {"targeting": [True, False], "element_id": ["a", "b"]},
        index=["g1", "g2"],
    )
    combined = pd.DataFrame({"targeting": [True, False]}, index=["g1", "g2"])

    restored = CONCAT_MODULE.preserve_source_guide_metadata(combined, source)

    assert str(restored["targeting"].dtype) == "bool"
    assert restored["element_id"].tolist() == ["a", "b"]
