
process prepare_all_guide_inference {

    cache 'lenient'

    input:
        path mudata
        path gtf_path

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input

    script:
        """
        create_pairs_to_test.py  --limit -1 ${mudata} ${gtf_path}
        prepare_inference.py pairs_to_test.csv ${mudata}
        """
}
