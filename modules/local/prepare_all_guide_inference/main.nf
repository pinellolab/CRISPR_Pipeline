
process prepare_all_guide_inference {

    cache 'lenient'

    input:
        path mudata
        path gtf_path

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input
        path "pairs_to_test.*", emit: pairs_to_test_file

    script:
        def pairs_format = ((params.INFERENCE_method == 'default' || params.INFERENCE_method.toString().contains('sceptre')) && params.INFERENCE_INTERMEDIATE_TABLE_FORMAT == 'parquet') ? 'tsv.gz' : params.INFERENCE_INTERMEDIATE_TABLE_FORMAT
        """
        prepare_inference.py \\
            --mudata_path ${mudata} \\
            --generate_pairs \\
            --input_gtf ${gtf_path} \\
            --limit -1 \\
            --pairs_output_format ${pairs_format}
        """
}
