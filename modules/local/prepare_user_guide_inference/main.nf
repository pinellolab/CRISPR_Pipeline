
process prepare_user_guide_inference {

    input:
        path mudata
        path user_inference

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input
        path "pairs_to_test.*", emit: pairs_to_test_file

    script:
        def pairs_format = ((params.INFERENCE_method == 'default' || params.INFERENCE_method.toString().contains('sceptre')) && params.INFERENCE_INTERMEDIATE_TABLE_FORMAT == 'parquet') ? 'tsv' : params.INFERENCE_INTERMEDIATE_TABLE_FORMAT
        """
        prepare_inference.py \\
            --guide_inference ${user_inference} \\
            --mudata_path ${mudata} \\
            --pairs_output_format ${pairs_format}
        """
}
