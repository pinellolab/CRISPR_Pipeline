"""Dependency-light regression tests for local/cis pair masking."""

from __future__ import annotations

import unittest

import anndata as ad
import jax.numpy as jnp
import numpy as np
from scipy import sparse

from perturbo.cli import _load_gene_by_element_mask
from perturbo.core import ControlFit, load_analysis_cells, subset_control_fit_genes
from perturbo.results import build_standard_element_effects_df


class LocalPairMaskTest(unittest.TestCase):
    def test_mask_is_aligned_to_requested_element_order(self) -> None:
        data = ad.AnnData(X=np.ones((2, 3), dtype=np.int32))
        data.varm["pairs"] = sparse.csr_matrix(
            np.asarray([[1, 0], [0, 1], [1, 0]], dtype=bool)
        )
        data.uns["pair_names"] = ["element_b", "element_a"]
        aligned = _load_gene_by_element_mask(
            data,
            modality_key=None,
            varm_key="pairs",
            names_uns_key="pair_names",
            element_names=["element_a", "element_b"],
        )
        np.testing.assert_array_equal(
            aligned.toarray(),
            np.asarray([[0, 1], [1, 0], [0, 1]], dtype=bool),
        )

    def test_analysis_loader_reads_only_selected_gene_columns(self) -> None:
        data = ad.AnnData(
            X=np.arange(12, dtype=np.int32).reshape(4, 3),
            obs={"perturbation": ["control", "a", "a", "control"]},
        )
        loaded = load_analysis_cells(
            data,
            perturbation_key="perturbation",
            selected_perturbations=["a"],
            selected_gene_indices=np.asarray([0, 2]),
            device="cpu",
        )
        self.assertEqual(loaded.counts.shape, (2, 2))
        np.testing.assert_array_equal(np.asarray(loaded.counts), np.asarray([[3, 5], [6, 8]]))

    def test_control_fit_gene_axes_are_subset(self) -> None:
        fit = ControlFit(
            beta_0=jnp.arange(4),
            theta=jnp.arange(4) + 10,
            noise_scale=jnp.ones(4),
            factor_loadings=jnp.arange(8).reshape(2, 1, 4),
            factor_scores=None,
            factor_center=jnp.arange(4),
            pca_loadings=jnp.arange(8).reshape(2, 1, 4),
            size_factors=jnp.ones((2, 1)),
            losses=jnp.ones(1),
            svi_result=None,
            covariate_coef=jnp.arange(8).reshape(2, 4),
        )
        selected = subset_control_fit_genes(fit, np.asarray([1, 3]))
        np.testing.assert_array_equal(np.asarray(selected.beta_0), np.asarray([1, 3]))
        np.testing.assert_array_equal(np.asarray(selected.factor_loadings), np.asarray([[[1, 3]], [[5, 7]]]))
        np.testing.assert_array_equal(np.asarray(selected.covariate_coef), np.asarray([[1, 3], [5, 7]]))

    def test_output_contains_only_masked_pairs(self) -> None:
        loc = np.arange(6, dtype=float).reshape(2, 3)
        scale = np.ones_like(loc)
        mask = sparse.csr_matrix(np.asarray([[1, 0, 1], [0, 1, 0]], dtype=bool))
        frame = build_standard_element_effects_df(
            method="perturbo",
            effect_loc=loc,
            effect_scale=scale,
            element_names=["a", "b"],
            gene_names=["g1", "g2", "g3"],
            tested_pairs_mask=mask,
        )
        self.assertEqual(
            set(zip(frame["element"], frame["gene"])),
            {("a", "g1"), ("a", "g3"), ("b", "g2")},
        )

    def test_output_without_mask_retains_logan_global_behavior(self) -> None:
        loc = np.arange(6, dtype=float).reshape(2, 3)
        scale = np.ones_like(loc)
        frame = build_standard_element_effects_df(
            method="perturbo",
            effect_loc=loc,
            effect_scale=scale,
            element_names=["a", "b"],
            gene_names=["g1", "g2", "g3"],
            tested_pairs_mask=None,
        )
        self.assertEqual(
            set(zip(frame["element"], frame["gene"])),
            {
                ("a", "g1"),
                ("a", "g2"),
                ("a", "g3"),
                ("b", "g1"),
                ("b", "g2"),
                ("b", "g3"),
            },
        )


if __name__ == "__main__":
    unittest.main()
