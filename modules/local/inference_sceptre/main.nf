
process inference_sceptre {

    input:
    path mudata_fp

    output:
    path "test_results.csv", emit: test_results

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
