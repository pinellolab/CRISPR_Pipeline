process INFERENCE_STUB {
    input:
    path concat_mudata
    path gtf_reference

    output:
    path 'inference_mudata.h5mu', emit: inference_mudata

    script:
    """
    cp ${concat_mudata} inference_mudata.h5mu
    """
}

workflow inference_pipeline {
    take:
    concat_mudata
    gtf_reference

    main:
    INFERENCE_STUB(concat_mudata, gtf_reference)

    emit:
    inference_mudata = INFERENCE_STUB.out.inference_mudata
}