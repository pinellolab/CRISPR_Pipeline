process postconcat_embedding_qc {
    tag "${stage}"
    cpus 4
    memory '32 GB'
    container { params.containers.base }
    publishDir { "${params.outdir}/postconcat_embedding_qc/${stage}" }, mode: params.publish_dir_mode, pattern: 'embedding_qc', overwrite: true

    input:
    path mudata_input
    val stage

    output:
    path 'postconcat_filtered.h5mu', emit: filtered_mudata
    path 'embedding_qc', emit: qc_dir

    script:
    """
    export MPLCONFIGDIR=\$PWD/.mpl NUMBA_CACHE_DIR=\$PWD/.numba
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    TAPSEQ_ARG='--no-tapseq-mode'
    if [ '${params.TAPSEQ_QC_MODE}' = 'true' ]; then TAPSEQ_ARG='--tapseq-mode'; fi
    LEIDEN_ARG='--no-enable-leiden'
    if [ '${params.QC_EMBEDDING_enable_leiden}' = 'true' ]; then LEIDEN_ARG='--enable-leiden'; fi
    python ${projectDir}/bin/postconcat_embedding_qc.py ${mudata_input} \
        --stage ${stage} --pct-mito ${params.QC_pct_mito} \
        --min-cells-fraction ${params.QC_min_cells_per_gene} \
        --n-pcs ${params.QC_EMBEDDING_n_pcs} --n-neighbors ${params.QC_EMBEDDING_n_neighbors} \
        --n-top-genes ${params.QC_EMBEDDING_n_top_genes} \
        --hvg-batch-key '${params.QC_EMBEDDING_hvg_batch_key}' \
        --batch-key '${params.QC_batch_col ?: 'batch'}' --seed ${params.QC_EMBEDDING_seed} \
        --max-dense-gb ${params.QC_EMBEDDING_max_dense_gb} \
        --reference '${params.REFERENCE_transcriptome}' \${TAPSEQ_ARG} \
        --cell-cycle '${params.QC_EMBEDDING_cell_cycle}' \
        --cell-cycle-genes ${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt \
        --min-cell-cycle-genes ${params.QC_EMBEDDING_min_cell_cycle_genes} \
        \${LEIDEN_ARG} --leiden-resolutions '${params.QC_EMBEDDING_leiden_resolutions}' \
        --leiden-diagnostic-resolution ${params.QC_EMBEDDING_leiden_diagnostic_resolution} \
        --leiden-n-iterations ${params.QC_EMBEDDING_leiden_n_iterations}
    """
}
