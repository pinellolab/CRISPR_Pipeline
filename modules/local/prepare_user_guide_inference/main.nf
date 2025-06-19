
process prepare_user_guide_inference {

    input:
        path mudata
        path user_inference

    output:
        path "mudata_inference_input.h5mu", emit: mudata_inference_input

    script:
        """
        prepare_inference.py ${user_inference} ${mudata} 
        """
}
