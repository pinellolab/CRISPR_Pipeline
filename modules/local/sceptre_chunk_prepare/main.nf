process sceptre_chunk_prepare {
    cache 'lenient'

    input:
    path mudata_fp

    output:
    path "chunks/*.h5mu", emit: mudata_chunks
    path "chunk_manifest.tsv", emit: chunk_manifest

    script:
    def force_chunk_flag = params.INFERENCE_SCEPTRE_FORCE_CHUNK ? '--force-chunk' : ''
    """
    chunk_mudata_sceptre.py \
        ${mudata_fp} \
        chunks \
        --chunk-by gene \
        --chunk-size ${params.INFERENCE_SCEPTRE_GENE_CHUNK_SIZE} \
        --chunk-mode ${params.INFERENCE_SCEPTRE_CHUNK_MODE} \
        --auto-threshold-entries ${params.INFERENCE_SCEPTRE_MAX_MATRIX_ENTRIES} \
        --keep-all-guides \
        --manifest chunk_manifest.tsv \
        ${force_chunk_flag}
    """
}
