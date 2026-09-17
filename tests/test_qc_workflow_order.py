from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_post_concat_mito_precedes_gene_support():
    script = (ROOT / "bin" / "concat_preprocessed_rna.py").read_text()
    assert script.index("mito_keep =") < script.index("required_cells = resolve_min_cells")


def test_clone_calling_is_after_aggregation_and_before_inference():
    workflow = (ROOT / "workflows" / "crispr_pipeline" / "main.nf").read_text()
    assert workflow.index("GuideAssignment = guide_assignment_pipeline") < workflow.index(
        "CloneRemoval = remove_clonal_cells"
    )
    assert workflow.index("CloneRemoval = remove_clonal_cells") < workflow.index(
        "Inference = inference_pipeline"
    )
    assert "CloneRemoval = remove_clonal_cells(GuideAssignment.concat_mudata)" in workflow


def test_scrublet_is_not_reapplied_after_mudata_aggregation():
    workflow = (ROOT / "workflows" / "crispr_pipeline" / "main.nf").read_text()
    assert "doublets_scrub(" not in workflow
    preprocessing = (ROOT / "subworkflows" / "local" / "preprocessing_pipeline" / "main.nf").read_text()
    assert "params.ENABLE_SCRUBLET" in preprocessing
    assert "params.SCRUBLET_assay_type" in preprocessing
