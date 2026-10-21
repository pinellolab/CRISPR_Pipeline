"""The non-targeting calibration input SCEPTRE receives.

Both regimes are built from one small screen: two targeting elements of two
guides each aimed at GENE1 and GENE2, four non-targeting guides that the
standard bucketing turns into two buckets, low-MOI assignments with a few
multi-guide cells, and the covariate columns a prepared MuData carries.
"""

import json
import pathlib
import sys

import anndata as ad
import h5py
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import prepare_ntc_calibration as prep
from intended_target_key_utils import _control_mask
from perturbo_v2_pipeline_adapter import _guide_control_mask

GENES = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6"]
GUIDES = ["g1", "g2", "g3", "g4", "nt1", "nt2", "nt3", "nt4"]
N_CELLS = 30


def _screen(seed=3):
    rng = np.random.default_rng(seed)
    cells = [f"cell{i:02d}" for i in range(N_CELLS)]
    # Rows: 0-9 carry one targeting guide each, 10-25 carry one NT guide each,
    # 26-27 carry a targeting and an NT guide, 28 carries two targeting guides,
    # 29 carries nothing.
    assignment = np.zeros((N_CELLS, len(GUIDES)), dtype=np.float32)
    for i in range(10):
        assignment[i, i % 4] = 1
    for i in range(10, 26):
        assignment[i, 4 + (i % 4)] = 1
    assignment[26, 0] = 1; assignment[26, 4] = 1
    assignment[27, 1] = 1; assignment[27, 5] = 1
    assignment[28, 2] = 1; assignment[28, 3] = 1
    guide_var = pd.DataFrame(
        {
            "guide_id": GUIDES,
            "targeting": [True] * 4 + [False] * 4,
            "type": ["targeting"] * 4 + ["non-targeting"] * 4,
            "intended_target_name": ["GENE1", "GENE1", "GENE2", "GENE2", "non-targeting", "non-targeting", "non-targeting", "non-targeting"],
            "intended_target_chr": ["chr1", "chr1", "chr2", "chr2", "", "", "", ""],
            "intended_target_start": [100.0, 100.0, 300.0, 300.0, np.nan, np.nan, np.nan, np.nan],
            "intended_target_end": [200.0, 200.0, 400.0, 400.0, np.nan, np.nan, np.nan, np.nan],
        },
        index=GUIDES,
    )
    counts = rng.poisson(5, size=(N_CELLS, len(GENES))).astype(np.float32) + 1
    gene_obs = pd.DataFrame(
        {
            "total_guide_umis": rng.integers(5, 50, N_CELLS).astype(float),
            "total_gene_umis": counts.sum(axis=1),
            "num_expressed_genes": (counts > 0).sum(axis=1).astype(float),
            "batch": np.where(np.arange(N_CELLS) % 2 == 0, "b1", "b2"),
        },
        index=cells,
    )
    for column in ("total_guide_umis", "total_gene_umis", "num_expressed_genes"):
        gene_obs[f"log_{column}"] = np.log(gene_obs[column])
    gene = ad.AnnData(X=sparse.csr_matrix(counts), obs=gene_obs, var=pd.DataFrame({"symbol": GENES}, index=GENES))
    guide = ad.AnnData(
        X=sparse.csr_matrix(assignment * 7), obs=pd.DataFrame(index=cells), var=guide_var,
    )
    guide.layers["guide_assignment"] = sparse.csr_matrix(assignment)
    guide.uns["moi"] = np.array(["low"], dtype=object)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    # What the prepared MuData carries: the design as top-level obs, the configured names in uns.
    mdata.obs = gene_obs[["total_guide_umis", "total_gene_umis", "num_expressed_genes", "batch"]].copy()
    mdata.uns["inference_continuous_covariates"] = ["log_total_guide_umis", "log_total_gene_umis", "log_num_expressed_genes"]
    mdata.uns["inference_batch_covariate"] = "batch"
    return mdata


def _pairs_mudata(mdata):
    pairs = pd.DataFrame(
        {
            "guide_id": ["g1", "g1", "g2", "g2", "g3", "g4"],
            "gene_id": ["GENE1", "GENE3", "GENE1", "GENE3", "GENE2", "GENE2"],
            "intended_target_key": [
                "coords::GENE1::chr1::100::200", "coords::GENE1::chr1::100::200",
                "coords::GENE1::chr1::100::200", "coords::GENE1::chr1::100::200",
                "coords::GENE2::chr2::300::400", "coords::GENE2::chr2::300::400",
            ],
            "pair_type": ["discovery"] * 6,
        }
    )
    out = mu.MuData({"gene": mdata["gene"].copy(), "guide": mdata["guide"].copy()})
    out.uns["pairs_to_test"] = pairs
    return out


def _stored_obs_columns(path):
    with h5py.File(path, "r") as handle:
        order = handle["obs"].attrs["column-order"]
    return [n.decode() if isinstance(n, bytes) else str(n) for n in order]


def _run(tmp_path, regime, **overrides):
    tmp_path.mkdir(parents=True, exist_ok=True)
    full_path = tmp_path / "full.h5mu"
    pairs_path = tmp_path / "pairs.h5mu"
    mdata = _screen()
    mdata.write(full_path)
    _pairs_mudata(mdata).write(pairs_path)
    out = {
        "output_mudata": tmp_path / f"{regime}_input.h5mu",
        "pairs": tmp_path / f"{regime}_pairs.tsv.gz",
        "guide_targets": tmp_path / f"{regime}_targets.tsv",
        "cell_counts_per_guide": tmp_path / f"{regime}_cells_guide.tsv",
        "cell_counts_per_element": tmp_path / f"{regime}_cells_element.tsv",
        "summary": tmp_path / f"{regime}_design.json",
    }
    argv = [
        "--full-mudata", str(full_path), "--pairs-mudata", str(pairs_path),
        "--sceptre-control-group", regime, "--analysis-moi", "low",
        "--max-local-pairs", str(overrides.get("max_local_pairs", 10000)),
        "--max-global-pairs", str(overrides.get("max_global_pairs", 10000)),
        "--seed", str(overrides.get("seed", 0)),
    ]
    out["direct_targets"] = tmp_path / f"{regime}_direct.tsv"
    for flag, path in out.items():
        argv += [f"--{flag.replace('_', '-')}", str(path)]
    assert prep.main(argv) == 0
    return out, json.loads(out["summary"].read_text())


def test_direct_target_pairs_name_every_targeting_guide_on_its_gene(tmp_path):
    out, _ = _run(tmp_path, "complement")
    direct = pd.read_csv(out["direct_targets"], sep="\t")
    assert set(map(tuple, direct[["guide_id", "gene_id"]].to_numpy())) == {
        ("g1", "GENE1"), ("g2", "GENE1"), ("g3", "GENE2"), ("g4", "GENE2")
    }


def test_complement_regime_keeps_every_cell_and_renames_the_buckets(tmp_path):
    out, summary = _run(tmp_path, "complement")
    assert summary["status"] == "ok"
    assert summary["n_ntc_buckets"] == 2
    result = mu.read_h5mu(out["output_mudata"])

    # Every cell, the pair genes only, the four NT guides as targets, no pool copy.
    assert result.n_obs == N_CELLS
    assert set(result["guide"].var["guide_id"]) == {"nt1", "nt2", "nt3", "nt4"}
    assert set(result["gene"].var_names) == set(pd.read_csv(out["pairs"], sep="\t")["gene_id"])
    var = result["guide"].var
    assert (var["type"] == "ntc_calibration").all()
    assert var["intended_target_name"].str.startswith("ntc_calibration|").all()
    assert var["intended_target_key"].str.startswith("name_only::ntc_calibration|").all()
    assert not _control_mask(var).any()
    assert not _guide_control_mask(var).any()
    # The covariates SCEPTRE reads (the stored top-level obs, as R reads it)
    # are exactly the source's design columns, not the modality columns
    # mu.read_h5mu pulls in under a prefix.
    assert _stored_obs_columns(out["output_mudata"]) == ["total_guide_umis", "total_gene_umis", "num_expressed_genes", "batch"]
    assert result.uns["inference_batch_covariate"] == "batch"
    pairs = pd.DataFrame(result.uns["pairs_to_test"])
    assert (pairs["pair_type"].astype(str) == "ntc_calibration").all()
    assert set(pairs.columns) >= {"guide_id", "gene_id", "intended_target_key", "pair_set", "ntc_bucket"}


def test_pair_sets_are_built_as_specified(tmp_path):
    out, summary = _run(tmp_path, "complement", max_global_pairs=6)
    pairs = pd.read_csv(out["pairs"], sep="\t")
    element_pairs = pairs.drop_duplicates(["ntc_bucket", "gene_id"])
    sets = element_pairs["pair_set"].astype(str)

    # local: the two targeting elements' gene lists, each on one random bucket.
    local = element_pairs[sets.str.contains("local")]
    assert set(local["gene_id"]) == {"GENE1", "GENE3", "GENE2"}
    assert summary["local_elements_used"] == 2
    # global: two buckets, cap 6 -> three sampled genes each.
    assert summary["global_genes_sampled"] == 3
    assert summary["global_pairs"] == 6
    # direct_target: both buckets against both direct-target genes.
    direct = element_pairs[sets.str.contains("direct_target")]
    assert set(direct["gene_id"]) == {"GENE1", "GENE2"}
    assert len(direct) == 4
    # Every guide of a bucket carries the bucket's pairs.
    per_bucket = pairs.groupby("ntc_bucket")["guide_id"].nunique()
    assert (per_bucket == 2).all()
    assert summary["union_guide_pairs"] == 2 * summary["union_element_pairs"]


def test_local_cap_keeps_whole_elements(tmp_path):
    out, summary = _run(tmp_path, "complement", max_local_pairs=2)
    # Elements hold 2 and 1 pairs; under a cap of 2 whole elements are kept in
    # random order while they fit, never a partial element.
    assert summary["local_pairs"] in (1, 2)
    assert summary["local_elements_used"] in (1, 2)
    pairs = pd.read_csv(out["pairs"], sep="\t")
    local = pairs[pairs["pair_set"].astype(str).str.contains("local")].drop_duplicates(["ntc_bucket", "gene_id"])
    assert len(local) == summary["local_pairs"]


def test_nt_cells_regime_duplicates_the_pool_with_a_synthetic_control(tmp_path):
    out, summary = _run(tmp_path, "nt_cells")
    assert summary["status"] == "ok"
    # 16 NT-only singleton cells + 2 NT-carrying multi-guide cells are kept; the
    # 16 NT-only cells are copied.
    assert summary["n_nt_only_cells"] == 16
    assert summary["n_pool_copies"] == 16
    assert summary["n_cells_kept"] == 18
    result = mu.read_h5mu(out["output_mudata"])
    assert result.n_obs == 18 + 16

    var = result["guide"].var
    assert prep.POOL_GUIDE_ID in set(var["guide_id"])
    pool = var.loc[var["guide_id"] == prep.POOL_GUIDE_ID].iloc[0]
    assert str(pool["intended_target_key"]).find("non-targeting") >= 0
    calibration = var.loc[var["guide_id"] != prep.POOL_GUIDE_ID]
    assert not _control_mask(calibration).any()
    assert _control_mask(var).sum() == 1

    assignment = result["guide"].layers["guide_assignment"].toarray()
    is_copy = pd.Index(result.obs_names).str.endswith(prep.POOL_CELL_SUFFIX)
    pool_col = list(var["guide_id"]).index(prep.POOL_GUIDE_ID)
    # Copies carry the pool gRNA and nothing else; originals never carry it.
    assert (assignment[is_copy, pool_col] == 1).all()
    assert (assignment[is_copy][:, [c for c in range(assignment.shape[1]) if c != pool_col]] == 0).all()
    assert (assignment[~is_copy, pool_col] == 0).all()
    # A copy is its original: same counts, same covariates.
    originals = [name[: -len(prep.POOL_CELL_SUFFIX)] for name in result.obs_names[is_copy]]
    gene = result["gene"]
    np.testing.assert_array_equal(gene[is_copy].X.toarray(), gene[originals].X.toarray())
    design = ["total_guide_umis", "total_gene_umis", "num_expressed_genes", "batch"]
    pd.testing.assert_frame_equal(
        result.obs.loc[is_copy, design].reset_index(drop=True), result.obs.loc[originals, design].reset_index(drop=True)
    )
    assert _stored_obs_columns(out["output_mudata"]) == ["total_guide_umis", "total_gene_umis", "num_expressed_genes", "batch"]
    ratios = summary["self_inclusion_ratio_by_bucket"]
    assert set(ratios) == {"1", "2"} and all(0 < r < 1 for r in ratios.values())


def test_same_seed_same_pairs(tmp_path):
    out_a, _ = _run(tmp_path / "a", "complement", seed=7)
    out_b, _ = _run(tmp_path / "b", "complement", seed=7)
    pd.testing.assert_frame_equal(
        pd.read_csv(out_a["pairs"], sep="\t"), pd.read_csv(out_b["pairs"], sep="\t")
    )


def test_screen_without_controls_is_skipped_cleanly(tmp_path):
    mdata = _screen()
    keep = mdata["guide"].var["targeting"].to_numpy()
    gene = mdata["gene"].copy()
    guide = mdata["guide"][:, keep].copy()
    targeting_only = mu.MuData({"gene": gene, "guide": guide})
    targeting_only.obs = mdata.obs.copy()
    full_path = tmp_path / "full.h5mu"
    targeting_only.write(full_path)
    pairs_path = tmp_path / "pairs.h5mu"
    _pairs_mudata(mdata).write(pairs_path)
    summary_path = tmp_path / "design.json"
    argv = [
        "--full-mudata", str(full_path), "--pairs-mudata", str(pairs_path),
        "--sceptre-control-group", "complement", "--summary", str(summary_path),
        "--output-mudata", str(tmp_path / "x.h5mu"), "--pairs", str(tmp_path / "x.tsv.gz"),
        "--guide-targets", str(tmp_path / "t.tsv"), "--direct-targets", str(tmp_path / "d.tsv"),
        "--cell-counts-per-guide", str(tmp_path / "cg.tsv"), "--cell-counts-per-element", str(tmp_path / "ce.tsv"),
    ]
    assert prep.main(argv) == 0
    summary = json.loads(summary_path.read_text())
    assert summary["status"] == "skipped"
    assert not (tmp_path / "x.h5mu").exists()
    assert (tmp_path / "cg.tsv").exists() and (tmp_path / "t.tsv").exists() and (tmp_path / "d.tsv").exists()


def test_cell_counts_report_roles_before_any_restriction(tmp_path):
    out, _ = _run(tmp_path, "nt_cells")
    per_guide = pd.read_csv(out["cell_counts_per_guide"], sep="\t")
    per_element = pd.read_csv(out["cell_counts_per_element"], sep="\t")
    assert set(per_guide["role"]) == {"targeting", "ntc"}
    # nt1 is carried by cells 10, 14, 18, 22 and the multi-guide cell 26.
    assert int(per_guide.loc[per_guide["guide_id"] == "nt1", "n_cells"].iloc[0]) == 5
    assert int(per_guide.loc[per_guide["guide_id"] == "nt1", "n_cells_singleton"].iloc[0]) == 4
    assert set(per_element["role"]) == {"targeting", "ntc"}
    assert per_element.loc[per_element["role"] == "ntc", "ntc_bucket"].notna().all()
