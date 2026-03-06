process sceptre_chunk_merge {
    cache 'lenient'

    input:
    path per_guide_outputs, stageAs: "guide?/*"
    path per_element_outputs, stageAs: "element?/*"
    path base_mudata
    path chunk_manifest
    val write_mudata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata, optional: true
    path "sceptre_per_element_output.tsv.gz", emit: per_element_output
    path "sceptre_per_guide_output.tsv.gz", emit: per_guide_output
    path "sceptre_merge_done.token", emit: completion_token

    script:
    def skipMudataArg = write_mudata ? '' : '--skip_mudata_write'
    """
    merge_sceptre_chunk_results.py \
        --per_guide_files ${per_guide_outputs.join(' ')} \
        --per_element_files ${per_element_outputs.join(' ')} \
        --base_mudata ${base_mudata} \
        --chunk_manifest ${chunk_manifest} \
        --output_mudata inference_mudata.h5mu \
        --output_per_guide sceptre_per_guide_output.tsv.gz \
        --output_per_element sceptre_per_element_output.tsv.gz \
        ${skipMudataArg}
    touch sceptre_merge_done.token
    """
}
