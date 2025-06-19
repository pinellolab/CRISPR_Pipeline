
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
    ${params.side}
    ${params.grna_integration_strategy}
    ${params.resampling_approximation}
    ${params.control_group}
    ${params.resampling_mechanism}
    ${params.formula_object}
    \$cov_string
    EOF

    inference_sceptre.R args.txt
    """
}
