
process mergedResults {
    cache 'lenient'
    publishDir './pipeline_outputs'

    input:
    path test_result
    path mudata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv", emit: per_element_output
    path "per_guide_output.tsv", emit: per_guide_output


    script:
        """
        export_output_multiple.py --sceptre_result ${test_result} --perturbo_mudata ${mudata}
        """

}
