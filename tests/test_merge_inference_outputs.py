import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import merge_cis_trans_results
import merge_method_results
import merge_sceptre_chunk_results


def _make_test_mudata():
    obs = pd.DataFrame(index=["cell1", "cell2"])
    gene = ad.AnnData(
        X=np.ones((2, 2), dtype=float),
        obs=obs.copy(),
        var=pd.DataFrame(index=["GENE1", "GENE2"]),
    )
    guide = ad.AnnData(
        X=np.ones((2, 2), dtype=float),
        obs=obs.copy(),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    return mu.MuData({"gene": gene, "guide": guide})


def _write_tsv(df, path):
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def test_merge_method_results_preserves_sceptre_se_and_adds_fdr(tmp_path, monkeypatch):
    base_mudata = tmp_path / "base.h5mu"
    _make_test_mudata().write(base_mudata)

    sceptre_guide = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "log2_fc": [1.0, -1.0],
            "p_value": [0.01, 0.2],
            "se_fold_change": [0.1, 0.2],
        }
    )
    sceptre_element = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "intended_target_name": ["target1", "target2"],
            "intended_target_chr": ["chr1", "chr2"],
            "intended_target_start": [100, 200],
            "intended_target_end": [150, 250],
            "log2_fc": [1.0, -1.0],
            "p_value": [0.01, 0.2],
            "se_fold_change": [0.1, 0.2],
        }
    )
    perturbo_guide = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "log2_fc": [0.5, -0.5],
            "p_value": [0.01, 0.2],
        }
    )
    perturbo_element = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "intended_target_name": ["target1", "target2"],
            "intended_target_chr": ["chr1", "chr2"],
            "intended_target_start": [100, 200],
            "intended_target_end": [150, 250],
            "log2_fc": [0.5, -0.5],
            "p_value": [0.01, 0.2],
        }
    )

    paths = {
        "sceptre_guide": tmp_path / "sceptre_guide.tsv.gz",
        "sceptre_element": tmp_path / "sceptre_element.tsv.gz",
        "perturbo_guide": tmp_path / "perturbo_guide.tsv.gz",
        "perturbo_element": tmp_path / "perturbo_element.tsv.gz",
    }
    _write_tsv(sceptre_guide, paths["sceptre_guide"])
    _write_tsv(sceptre_element, paths["sceptre_element"])
    _write_tsv(perturbo_guide, paths["perturbo_guide"])
    _write_tsv(perturbo_element, paths["perturbo_element"])

    monkeypatch.chdir(tmp_path)
    merge_method_results.merge_method_results(
        str(paths["sceptre_guide"]),
        str(paths["sceptre_element"]),
        str(paths["perturbo_guide"]),
        str(paths["perturbo_element"]),
        str(base_mudata),
    )

    guide_out = pd.read_csv(tmp_path / "per_guide_output.tsv.gz", sep="\t")
    element_out = pd.read_csv(tmp_path / "per_element_output.tsv.gz", sep="\t")

    for observed in (guide_out, element_out):
        assert "sceptre_q_value" in observed.columns
        assert "sceptre_fc_se" in observed.columns
        assert "perturbo_fdr_log10_p_value" in observed.columns

    assert np.isclose(guide_out.loc[0, "sceptre_q_value"], 0.02)
    assert np.isclose(guide_out.loc[0, "sceptre_fc_se"], 0.1)
    assert np.isclose(
        guide_out.loc[0, "perturbo_fdr_log10_p_value"], -np.log10(0.02)
    )

    mdata = mu.read_h5mu(tmp_path / "inference_mudata.h5mu")
    stored = pd.DataFrame(mdata.uns["per_element_results"])
    assert "sceptre_q_value" in stored.columns
    assert "perturbo_fdr_log10_p_value" in stored.columns


def test_merge_sceptre_chunk_results_adds_global_q_values(tmp_path):
    base_mudata = tmp_path / "base.h5mu"
    _make_test_mudata().write(base_mudata)

    guide_chunk = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "log2_fc": [1.0, -1.0],
            "p_value": [0.01, 0.2],
        }
    )
    element_chunk = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "intended_target_name": ["target1", "target2"],
            "log2_fc": [1.0, -1.0],
            "p_value": [0.01, 0.2],
        }
    )
    guide_path = tmp_path / "guide.tsv.gz"
    element_path = tmp_path / "element.tsv.gz"
    manifest_path = tmp_path / "manifest.tsv"
    _write_tsv(guide_chunk, guide_path)
    _write_tsv(element_chunk, element_path)
    manifest_path.write_text("chunk_id\n0\n")

    merge_sceptre_chunk_results.merge_sceptre_chunk_results(
        per_guide_files=[str(guide_path)],
        per_element_files=[str(element_path)],
        base_mudata=str(base_mudata),
        chunk_manifest=str(manifest_path),
        output_mudata=str(tmp_path / "sceptre_mudata.h5mu"),
        output_per_guide=str(tmp_path / "sceptre_per_guide.tsv.gz"),
        output_per_element=str(tmp_path / "sceptre_per_element.tsv.gz"),
    )

    observed = pd.read_csv(tmp_path / "sceptre_per_guide.tsv.gz", sep="\t")
    assert "q_value" in observed.columns
    assert np.isclose(observed.loc[0, "q_value"], 0.02)


def test_merge_cis_trans_results_writes_fdr_columns_to_outputs_and_mudata(
    tmp_path, monkeypatch
):
    base_mudata = tmp_path / "base.h5mu"
    _make_test_mudata().write(base_mudata)

    cis_guide = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "sceptre_log2_fc": [1.0, -1.0],
            "sceptre_p_value": [0.01, 0.2],
            "sceptre_q_value": [0.02, 0.2],
            "sceptre_fc_se": [0.1, 0.2],
            "perturbo_log2_fc": [0.5, -0.5],
            "perturbo_p_value": [0.01, 0.2],
        }
    )
    cis_element = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "intended_target_name": ["target1", "target2"],
            "intended_target_chr": ["chr1", "chr2"],
            "intended_target_start": [100, 200],
            "intended_target_end": [150, 250],
            "sceptre_log2_fc": [1.0, -1.0],
            "sceptre_p_value": [0.01, 0.2],
            "sceptre_q_value": [0.02, 0.2],
            "sceptre_fc_se": [0.1, 0.2],
            "perturbo_log2_fc": [0.5, -0.5],
            "perturbo_p_value": [0.01, 0.2],
        }
    )
    trans_guide = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "log2_fc": [0.4, -0.4],
            "p_value": [0.05, 0.5],
        }
    )
    trans_element = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "intended_target_name": ["target1", "target2"],
            "intended_target_chr": ["chr1", "chr2"],
            "intended_target_start": [100, 200],
            "intended_target_end": [150, 250],
            "log2_fc": [0.4, -0.4],
            "p_value": [0.05, 0.5],
        }
    )

    paths = {
        "cis_guide": tmp_path / "cis_guide.tsv.gz",
        "cis_element": tmp_path / "cis_element.tsv.gz",
        "trans_guide": tmp_path / "trans_guide.tsv.gz",
        "trans_element": tmp_path / "trans_element.tsv.gz",
    }
    _write_tsv(cis_guide, paths["cis_guide"])
    _write_tsv(cis_element, paths["cis_element"])
    _write_tsv(trans_guide, paths["trans_guide"])
    _write_tsv(trans_element, paths["trans_element"])

    monkeypatch.chdir(tmp_path)
    merge_cis_trans_results.merge_cis_trans_results(
        str(paths["cis_guide"]),
        str(paths["cis_element"]),
        str(paths["trans_guide"]),
        str(paths["trans_element"]),
        str(base_mudata),
        str(tmp_path / "inference_mudata.h5mu"),
    )

    cis_observed = pd.read_csv(tmp_path / "cis_per_element_output.tsv.gz", sep="\t")
    trans_observed = pd.read_csv(tmp_path / "trans_per_element_output.tsv.gz", sep="\t")

    assert "sceptre_q_value" in cis_observed.columns
    assert "sceptre_fc_se" in cis_observed.columns
    assert "perturbo_fdr_log10_p_value" in cis_observed.columns
    assert "perturbo_fdr_log10_p_value" in trans_observed.columns
    assert np.isclose(
        trans_observed.loc[0, "perturbo_fdr_log10_p_value"], -np.log10(0.1)
    )

    mdata = mu.read_h5mu(tmp_path / "inference_mudata.h5mu")
    stored_cis = pd.DataFrame(mdata.uns["cis_per_element_results"])
    stored_trans = pd.DataFrame(mdata.uns["trans_per_element_results"])
    assert "sceptre_q_value" in stored_cis.columns
    assert "sceptre_fc_se" in stored_cis.columns
    assert "perturbo_fdr_log10_p_value" in stored_cis.columns
    assert "perturbo_fdr_log10_p_value" in stored_trans.columns
