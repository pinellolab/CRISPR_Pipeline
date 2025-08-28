
process inference_sceptre {

    input:
    path mudata_fp

    output:
    path "sceptre_inference_mudata.h5mu", emit: inference_mudata
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
    EOF

    inference_sceptre.R args.txt
    """
}
