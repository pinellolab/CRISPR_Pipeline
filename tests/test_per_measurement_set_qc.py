import importlib.util
from argparse import Namespace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


SCRIPT = Path(__file__).parents[1] / "bin" / "preprocess_adata.py"
SPEC = importlib.util.spec_from_file_location("preprocess_adata", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONCAT_SCRIPT = Path(__file__).parents[1] / "bin" / "concat_preprocessed_rna.py"
CONCAT_SPEC = importlib.util.spec_from_file_location("concat_preprocessed_rna", CONCAT_SCRIPT)
CONCAT_MODULE = importlib.util.module_from_spec(CONCAT_SPEC)
CONCAT_SPEC.loader.exec_module(CONCAT_MODULE)


def test_mad_limits_are_metric_specific():
    median, mad, lower, upper = MODULE.mad_limits([1, 2, 3, 4, 100], 3)
    assert median == 3
    assert mad == 1
    assert lower == 0
    assert upper == 6

    _, _, lower, upper = MODULE.mad_limits([1, 2, 3, 4, 100], 3, upper_only=True)
    assert np.isneginf(lower)
    assert upper == 6

    _, _, lower, upper = MODULE.mad_limits([1, 2, 3], 0)
    assert np.isneginf(lower)
    assert np.isposinf(upper)


def test_one_measurement_set_is_filtered_and_audited(tmp_path, monkeypatch):
    mapping = tmp_path / "B1_ks_transcripts_out"
    counts = mapping / "counts_unfiltered"
    counts.mkdir(parents=True)
    matrix = sparse.csr_matrix(
        np.array(
            [
                [10, 3, 5],
                [8, 2, 3],
                [0, 1, 0],
                [4, 2, 1],
            ],
            dtype=np.int32,
        )
    )
    raw = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=["AAAA", "AAAC", "AAAG", "AAAT"]),
        var=pd.DataFrame(index=["ENSG1.1", "ENSG2.1", "ENSG3.1"]),
    )
    raw.write_h5ad(counts / "adata.h5ad")
    (counts / "cells_x_genes.genes.names.txt").write_text("MT-CO1\nRPLP0\nGENE1\n")
    covariates = tmp_path / "parse_covariate.csv"
    pd.DataFrame(
        {"batch": ["B1"], "concat_batch": ["sample_0"], "barcode_key": ["B1"]}
    ).to_csv(covariates, index=False)

    monkeypatch.chdir(tmp_path)
    MODULE.main(
        Namespace(
            mapping_dir=str(mapping),
            covariates=str(covariates),
            qc_dir=tmp_path / "B1_qc",
            min_counts=2,
            reference="human",
            barcode_filter="none",
            mad_total_counts=0,
            mad_n_genes=0,
            enable_scrublet=False,
            scrublet_expected_doublet_rate=0.08,
            bc_replacement=False,
            use_multimapping=False,
        )
    )

    filtered = ad.read_h5ad(tmp_path / "B1_filtered.h5ad")
    assert filtered.n_obs == 3
    assert filtered.obs_names.tolist() == ["AAAA_B1", "AAAC_B1", "AAAT_B1"]
    assert filtered.obs["batch"].astype(str).unique().tolist() == ["B1"]

    audit = pd.read_csv(tmp_path / "B1_qc" / "measurement_set_qc_B1.tsv", sep="\t")
    assert audit.loc[0, "measurement_set"] == "B1"
    assert audit.loc[0, "input_barcodes"] == 4
    assert audit.loc[0, "retained_cells"] == 3
    assert (tmp_path / "B1_qc" / "knee_plot_scRNA_B1.png").exists()
    assert (tmp_path / "B1_qc" / "qc_distributions_scRNA_B1.png").exists()
    assert (tmp_path / "B1_qc" / "rna_qc_filter_flow_B1.png").exists()
    assert (tmp_path / "B1_qc" / "rna_qc_filter_steps_B1.png").exists()
    flow = pd.read_csv(tmp_path / "B1_qc" / "rna_qc_filter_flow_B1.tsv", sep="\t")
    assert flow["filter_label"].tolist() == [
        "QC_barcode_filter = none",
        "QC_min_counts_per_cell",
        "QC_MAD_total_counts",
        "QC_MAD_n_genes",
        "Scrublet doublet removal",
    ]
    count_step = flow.loc[flow["filter_label"] == "QC_min_counts_per_cell"].iloc[0]
    assert count_step["cells_before"] == 4
    assert count_step["cells_after"] == 3
    assert count_step["cells_removed"] == 1


def test_concatenated_gene_metrics_are_recomputed():
    matrix = sparse.csr_matrix(np.array([[1, 0], [2, 3]], dtype=np.int32))
    combined = ad.AnnData(X=matrix, var=pd.DataFrame({"symbol": ["A", "B"]}, index=["g1", "g2"]))

    detected = CONCAT_MODULE.recompute_gene_metrics(combined)

    assert detected.tolist() == [2, 1]
    assert combined.var["total_counts"].tolist() == [3, 3]
    assert combined.var["pct_dropout_by_counts"].tolist() == [0, 50]


def test_post_concat_mito_then_fractional_gene_filter(tmp_path, monkeypatch):
    var = pd.DataFrame(index=["g1", "g2"])
    first = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [2, 1]], dtype=np.int32),
        obs=pd.DataFrame({"pct_counts_mt": [5.0, 30.0]}, index=["c1", "c2"]),
        var=var.copy(),
    )
    second = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [2, 1]], dtype=np.int32),
        obs=pd.DataFrame({"pct_counts_mt": [6.0, 7.0]}, index=["c3", "c4"]),
        var=var.copy(),
    )
    first.write_h5ad(tmp_path / "a.h5ad")
    second.write_h5ad(tmp_path / "b.h5ad")
    monkeypatch.setattr(
        "sys.argv",
        [
            "concat_preprocessed_rna.py",
            str(tmp_path / "a.h5ad"),
            str(tmp_path / "b.h5ad"),
            "--output", str(tmp_path / "filtered.h5ad"),
            "--qc-dir", str(tmp_path / "post_qc"),
            "--pct-mito", "15",
            "--min-cells-fraction", "0.5",
        ],
    )
    CONCAT_MODULE.main()
    filtered = ad.read_h5ad(tmp_path / "filtered.h5ad")
    assert filtered.obs_names.tolist() == ["c1", "c3", "c4"]
    assert filtered.var_names.tolist() == ["g1"]
    flow = pd.read_csv(tmp_path / "post_qc" / "post_concat_qc_filter_flow.tsv", sep="\t")
    assert flow["filter"].tolist() == ["QC_pct_mito", "QC_min_cells_per_gene"]
    assert flow.loc[0, ["before", "after"]].tolist() == [4, 3]
    assert flow.loc[1, ["before", "after"]].tolist() == [2, 1]


def test_covariate_file_is_broadcast_to_every_measurement_set():
    workflow = (
        Path(__file__).parents[1] / "subworkflows" / "local" / "preprocessing_pipeline" / "main.nf"
    ).read_text()
    assert "parsed_covariate_value = parsed_covariate_file.first()" in workflow
    assert "PreprocessAnnData(\n        trans_out_dir,\n        parsed_covariate_value," in workflow
