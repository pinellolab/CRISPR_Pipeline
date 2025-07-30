
process inference_sceptre {

    input:
    path mudata_fp
    path cov_string

    output:
    path "test_results.csv", emit: test_results

    script:
    """
    cov_string=\$(cat $cov_string)

    cat <<EOF > args.txt
    ${mudata_fp}
    ${params.INFERENCE_SCEPTRE_side}
    ${params.INFERENCE_SCEPTRE_grna_integration_strategy}
    ${params.INFERENCE_SCEPTRE_resampling_approximation}
    ${params.INFERENCE_SCEPTRE_control_group}
    ${params.INFERENCE_SCEPTRE_resampling_mechanism}
    default
    \$cov_string
    EOF

    inference_sceptre.R args.txt
    """
}
