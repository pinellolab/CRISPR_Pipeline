process concat_preprocessed_rna {
    cache 'lenient'

    input:
    path filtered_measurement_sets
    val tapseq_qc_mode

    output:
    path "filtered_anndata.h5ad", emit: filtered_anndata_rna

    script:
    def inputs = filtered_measurement_sets.collect { it.toString() }.sort().join(' ')
    def tapseqArg = tapseq_qc_mode ? '--tapseq-mode' : ''
    """
    concat_preprocessed_rna.py ${inputs} --output filtered_anndata.h5ad ${tapseqArg}
    """
}
