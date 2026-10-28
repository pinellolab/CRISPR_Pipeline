process balanced_normalization {
    cache 'lenient'
    cpus 4
    memory '64 GB'
    time '24h'
    publishDir "${params.outdir}/postconcat_embedding_qc/balanced", mode: 'copy'
    input:
    path mudata
    path cell_cycle_genes
    output:
    path 'embedding_qc', emit: qc_dir
    script:
    """
    export MPLCONFIGDIR=\$PWD/.mplconfig
    python ${projectDir}/bin/balanced_normalization.py --input ${mudata} \\
        --batch-key '${params.QC_batch_col ?: 'batch'}' --seed ${params.BALANCED_seed} \\
        --cell-cycle-genes ${cell_cycle_genes}
    """
}
