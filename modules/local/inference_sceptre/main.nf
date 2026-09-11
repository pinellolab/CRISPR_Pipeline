
process inference_sceptre {

    input:
    path mudata_fp

    output:
    // Optional: nothing downstream reads a chunk's MuData (sceptre_chunk_merge takes
    // the tables), and on a screen-scale input each one is gigabytes.
    path "inference_mudata.h5mu", optional: true, emit: inference_mudata
    path "sceptre_per_element_output.tsv.gz", emit: per_element_output
    path "sceptre_per_guide_output.tsv.gz", emit: per_guide_output

    script:
    """
    cat <<EOF > args.txt
    ${mudata_fp}
    ${params.INFERENCE_SCEPTRE_side}
    ${params.INFERENCE_SCEPTRE_grna_integration_strategy}
    ${params.INFERENCE_SCEPTRE_resampling_approximation}
    ${params.INFERENCE_SCEPTRE_control_group}
    ${params.INFERENCE_SCEPTRE_resampling_mechanism}
    ${task.cpus}
    EOF

    export SCEPTRE_WRITE_MUDATA=${params.INFERENCE_SCEPTRE_WRITE_MUDATA}
    inference_sceptre.R args.txt
    mv per_element_output.tsv.gz sceptre_per_element_output.tsv.gz
    mv per_guide_output.tsv.gz sceptre_per_guide_output.tsv.gz
    """
}
