
process prepare_user_guide_inference {

    input:
        path mudata
        path user_inference

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input

    script:
        """
        prepare_inference.py ${user_inference} ${mudata} \
          --cov-mito ${params.INFERENCE_COV_MITOCHONDRIAL} --cov-batch ${params.INFERENCE_COV_BATCH} \
          --cov-ribo ${params.INFERENCE_COV_RIBOSOMAL} --cov-detected-genes ${params.INFERENCE_COV_DETECTED_GENES} \
          --cov-guide-count ${params.INFERENCE_COV_GUIDE_COUNT} --cov-guide-umis ${params.INFERENCE_COV_GUIDE_UMIS} \
          --cov-custom '${params.INFERENCE_COV_CUSTOM}'
        """
}
