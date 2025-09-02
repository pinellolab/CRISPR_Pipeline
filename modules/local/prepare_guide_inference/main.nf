
process prepare_guide_inference {

    cache 'lenient'

    input:
        path mudata
        path gtf_path
        val limit
        val subset_for_cis

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input

    script:
        def cis_flag = subset_for_cis ? "--subset_for_cis" : ""
        """
        create_pairs_to_test.py  --limit $limit ${mudata} ${gtf_path}
        prepare_inference.py pairs_to_test.csv ${mudata} ${cis_flag}
        """
}
