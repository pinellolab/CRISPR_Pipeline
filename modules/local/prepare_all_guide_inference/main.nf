
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
        prepare_inference.py pairs_to_test.csv ${mudata} \
          --cov-mito ${params.INFERENCE_COV_MITOCHONDRIAL} --cov-batch ${params.INFERENCE_COV_BATCH} \
          --cov-ribo ${params.INFERENCE_COV_RIBOSOMAL} --cov-detected-genes ${params.INFERENCE_COV_DETECTED_GENES} \
          --cov-guide-count ${params.INFERENCE_COV_GUIDE_COUNT} --cov-guide-umis ${params.INFERENCE_COV_GUIDE_UMIS} \
          --cov-custom '${params.INFERENCE_COV_CUSTOM}'
        """
}
