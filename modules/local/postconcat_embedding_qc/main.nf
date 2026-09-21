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
    python ${projectDir}/bin/postconcat_embedding_qc.py ${mudata_input} \
        --stage ${stage} --pct-mito ${params.QC_pct_mito} \
        --min-cells-fraction ${params.QC_min_cells_per_gene} \
        --n-pcs ${params.QC_EMBEDDING_n_pcs} --n-neighbors ${params.QC_EMBEDDING_n_neighbors} \
        --n-top-genes ${params.QC_EMBEDDING_n_top_genes} \
        --hvg-batch-key '${params.QC_EMBEDDING_hvg_batch_key}' \
        --batch-key '${params.QC_batch_col ?: 'batch'}' --seed ${params.QC_EMBEDDING_seed} \
        --max-dense-gb ${params.QC_EMBEDDING_max_dense_gb}
    """
}
