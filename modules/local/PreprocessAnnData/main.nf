
process PreprocessAnnData {

    cache 'lenient'
    tag { mapping_dir.getBaseName().replaceFirst(/_ks_transcripts_out$/, '') }

    input:
    path mapping_dir
    path parsed_covariates
    val min_genes
    val min_counts
    val pct_mito
    val reference
    val barcode_filter
    val mad_total_counts
    val mad_n_genes
    val mad_pct_mito

    output:
    path "*_filtered.h5ad", emit: filtered_measurement_set
    path "*_qc", emit: measurement_set_qc

    script:
        def batch = mapping_dir.getBaseName().replaceFirst(/_ks_transcripts_out$/, '')
        def safeBatch = batch.replaceAll(/[^A-Za-z0-9_.-]+/, '_')
        def bcArg = params.replace_barcodes ? '--bc-replacement' : ''
        def mmArg = params.use_multimapping ? '--use-multimapping' : ''
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}
        preprocess_adata.py \
            ${mapping_dir} ${parsed_covariates} \
            --qc-dir ${safeBatch}_qc \
            --min-genes ${min_genes} \
            --min-counts ${min_counts} \
            --pct-mito ${pct_mito} \
            --reference ${reference} \
            --barcode-filter ${barcode_filter} \
            --mad-total-counts ${mad_total_counts} \
            --mad-n-genes ${mad_n_genes} \
            --mad-pct-mito ${mad_pct_mito} \
            ${bcArg} ${mmArg}
        """
}
