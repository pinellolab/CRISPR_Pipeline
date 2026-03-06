import os
import sys
import unittest

import pandas as pd

BIN_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "bin")
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)

from intended_target_key_utils import annotate_intended_target_groups, enrich_pairs_with_target_metadata


class IntendedTargetKeyUtilsTests(unittest.TestCase):
    def test_distinct_coordinate_keys_for_same_name(self):
        guide_var = pd.DataFrame(
            {
                "guide_id": ["g1", "g2"],
                "targeting": [True, True],
                "type": ["targeting", "targeting"],
                "intended_target_name": ["ENSG1", "ENSG1"],
                "intended_target_chr": ["chr1", "chr1"],
                "intended_target_start": [100, 200],
                "intended_target_end": [110, 210],
            }
        )

        out = annotate_intended_target_groups(guide_var)
        self.assertEqual(out["intended_target_key"].nunique(), 2)

    def test_non_targeting_are_bucketed_and_never_exact_non_targeting(self):
        guide_var = pd.DataFrame(
            {
                "guide_id": ["g1", "g2", "g3", "g4", "nt1", "nt2"],
                "targeting": [True, True, True, True, False, False],
                "type": ["targeting", "targeting", "targeting", "targeting", "non-targeting", "non-targeting"],
                "intended_target_name": ["E1", "E1", "E2", "E3", "non-targeting", "non-targeting"],
                "intended_target_chr": ["chr1", "chr1", "chr2", "chr3", pd.NA, pd.NA],
                "intended_target_start": [10, 10, 20, 30, pd.NA, pd.NA],
                "intended_target_end": [11, 11, 21, 31, pd.NA, pd.NA],
            }
        )

        out = annotate_intended_target_groups(guide_var)
        self.assertFalse((out["intended_target_name"] == "non-targeting").any())
        control_rows = out[out["guide_id"].isin(["nt1", "nt2"])]
        self.assertTrue(control_rows["intended_target_name"].str.startswith("non-targeting|").all())

    def test_missing_coords_falls_back_to_name_only_key(self):
        guide_var = pd.DataFrame(
            {
                "guide_id": ["g1"],
                "targeting": [True],
                "type": ["targeting"],
                "intended_target_name": ["ENSGX"],
                "intended_target_chr": [pd.NA],
                "intended_target_start": [pd.NA],
                "intended_target_end": [pd.NA],
            }
        )

        out = annotate_intended_target_groups(guide_var)
        self.assertTrue(out.loc[0, "intended_target_key"].startswith("name_only::ENSGX"))

    def test_enrich_pairs_autofills_coordinates_and_key(self):
        guide_var = pd.DataFrame(
            {
                "guide_id": ["g1"],
                "targeting": [True],
                "type": ["targeting"],
                "intended_target_name": ["ENSG1"],
                "intended_target_chr": ["chr1"],
                "intended_target_start": [100],
                "intended_target_end": [110],
            }
        )
        guide_var = annotate_intended_target_groups(guide_var)

        pairs = pd.DataFrame(
            {
                "guide_id": ["g1"],
                "gene_name": ["ENSG1"],
                "pair_type": ["discovery"],
            }
        )

        out = enrich_pairs_with_target_metadata(pairs, guide_var)
        self.assertEqual(out.loc[0, "intended_target_chr"], "chr1")
        self.assertEqual(int(out.loc[0, "intended_target_start"]), 100)
        self.assertEqual(int(out.loc[0, "intended_target_end"]), 110)
        self.assertTrue(isinstance(out.loc[0, "intended_target_key"], str))
        self.assertGreater(len(out.loc[0, "intended_target_key"]), 0)


if __name__ == "__main__":
    unittest.main()
