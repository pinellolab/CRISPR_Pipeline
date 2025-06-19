process inference_mudata {
    cache 'lenient'
    publishDir './pipeline_outputs'

    input:
    path test_result
    path mudata
    val inference_method

    output:
        path "inference_mudata.h5mu", emit: inference_mudata
        path "per_element_output.tsv", emit: per_element_output
        path "per_guide_output.tsv", emit: per_guide_output

    script:
        """
        add_guide_inference.py --test_results_csv ${test_result} --mudata ${mudata}
        export_output_single.py --mudata inference_mudata.h5mu --inference_method ${inference_method}
        """

}
