
process inference_sceptre {

    input:
    path mudata_fp
    path pairs_to_test_fp

    output:
    path "sceptre_per_element_output.tsv.gz", emit: per_element_output
    path "sceptre_per_guide_output.tsv.gz", emit: per_guide_output

    script:
    def pairsArg = pairs_to_test_fp.getName() == 'no_pairs_to_test.tsv' ? 'NONE' : pairs_to_test_fp
    """
    cat <<EOF > args.txt
    ${mudata_fp}
    ${params.INFERENCE_SCEPTRE_side}
    ${params.INFERENCE_SCEPTRE_grna_integration_strategy}
    ${params.INFERENCE_SCEPTRE_resampling_approximation}
    ${params.INFERENCE_SCEPTRE_control_group}
    ${params.INFERENCE_SCEPTRE_resampling_mechanism}
    ${task.cpus}
    ${pairsArg}
    EOF

    inference_sceptre.R args.txt
    mv per_element_output.tsv.gz sceptre_per_element_output.tsv.gz
    mv per_guide_output.tsv.gz sceptre_per_guide_output.tsv.gz
    """
}
