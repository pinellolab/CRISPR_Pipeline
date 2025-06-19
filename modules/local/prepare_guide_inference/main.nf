
process prepare_guide_inference {

    cache 'lenient'

    input:
        path mudata
        path gtf_path
        val limit

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input

    script:
        """
        create_pairs_to_test.py  --limit $limit ${mudata} ${gtf_path}
        prepare_inference.py pairs_to_test.csv ${mudata}
        """
}
