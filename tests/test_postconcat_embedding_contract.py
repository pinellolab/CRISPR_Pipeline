from pathlib import Path

ROOT = Path(__file__).parents[1]

def test_pre_and_post_clone_use_independent_raw_sources():
    workflow = (ROOT/'workflows/crispr_pipeline/main.nf').read_text()
    assert "embedding_before_clone(qualified_mudata, 'before_clone', file(" in workflow
    assert 'remove_clonal_cells(qualified_mudata)' in workflow
    assert "embedding_after_clone(mudata_before_hto, 'after_clone', file(" in workflow
    assert 'remove_clonal_cells(BeforeEmbedding.filtered_mudata)' not in workflow


def test_cell_cycle_markers_are_staged_as_a_process_input():
    module = (ROOT/"modules/local/postconcat_embedding_qc/main.nf").read_text()
    workflow = (ROOT/"workflows/crispr_pipeline/main.nf").read_text()
    assert "path cell_cycle_genes" in module
    assert "--cell-cycle-genes ${cell_cycle_genes}" in module
    assert "${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt" in workflow


def test_counts_not_normalized_in_place():
    script = (ROOT/'bin/postconcat_embedding_qc.py').read_text()
    assert 'sc.pp.normalize_total(temp, target_sum=None)' in script
    assert "sp.csr_matrix(raw.X, dtype=np.float32).copy()" in script
    assert 'filtered.write_h5mu(args.output)' in script
    assert 'temp.write' not in script
    assert 'faceted_measurement_pca(temp, args)' in script
    assert "result['cell_cycle'] = score_cell_cycle(temp, args)" in script
    assert "result['leiden'] = run_leiden_sweep(temp, args)" in script
    assert "temp.obs[[args.batch_key] + keys].to_csv" in script

def test_concat_defers_filters():
    script = (ROOT/'bin/concat_preprocessed_rna.py').read_text()
    assert script.index('if args.defer_global_qc:') < script.index('mito_keep =')
