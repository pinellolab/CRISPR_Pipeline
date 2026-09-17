process concat_preprocessed_rna {
    cache 'lenient'

    input:
    path filtered_measurement_sets
    val pct_mito
    val min_cells_fraction

    output:
    path "filtered_anndata.h5ad", emit: filtered_anndata_rna
    path "post_concat_qc", emit: post_concat_qc

    script:
    def inputs = filtered_measurement_sets.collect { it.toString() }.sort().join(' ')
    """
    concat_preprocessed_rna.py ${inputs} \
        --output filtered_anndata.h5ad \
        --qc-dir post_concat_qc \
        --pct-mito ${pct_mito} \
        --min-cells-fraction ${min_cells_fraction}
    """
}
