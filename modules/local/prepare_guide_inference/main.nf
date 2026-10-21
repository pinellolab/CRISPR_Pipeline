
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
        prepare_inference.py pairs_to_test.csv ${mudata} ${cis_flag} \
          --cov-mito ${params.INFERENCE_COV_MITOCHONDRIAL} --cov-batch ${params.INFERENCE_COV_BATCH} \
          --cov-ribo ${params.INFERENCE_COV_RIBOSOMAL} --cov-detected-genes ${params.INFERENCE_COV_DETECTED_GENES} \
          --cov-guide-count ${params.INFERENCE_COV_GUIDE_COUNT} --cov-guide-umis ${params.INFERENCE_COV_GUIDE_UMIS} \
          --cov-custom '${params.INFERENCE_COV_CUSTOM}'
        """
}
