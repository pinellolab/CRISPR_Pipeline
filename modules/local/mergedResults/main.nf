
process mergedResults {
    cache 'lenient'
    publishDir './pipeline_outputs', mode: 'copy', overwrite: true

    input:
    path test_result
    path mudata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv", emit: per_element_output
    path "per_guide_output.tsv", emit: per_guide_output


    script:
        """
        # Rename input to avoid conflicts
        [[ -e ${mudata} ]] && mv ${mudata} input_mudata.h5mu

        export_output_multiple.py --sceptre_result ${test_result} --perturbo_mudata input_mudata.h5mu
        """

}
