
process inference_perturbo {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, enabled: { params.INFERENCE_method == 'perturbo' }

    input:
    path mudata
    val inference_method

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_local_analysis_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_local_analysis_per_guide_output.tsv.gz", emit: per_guide_output
    path "perturbo_v2_outputs", optional: true, emit: perturbo_v2_outputs

    script:
        def save_model_params_arg = params.INFERENCE_PERTURBO_SAVE_MODEL_PARAMS ? '--save-model-params' : '--no-save-model-params'
        def parallel_fits_arg = params.INFERENCE_PERTURBO_LOCAL_PARALLEL_FITS ? '--parallel-fits' : '--no-parallel-fits'
        """
        perturbo_v2_pipeline_adapter.py \\
            --input ${mudata} \\
            --per-element-output perturbo_local_analysis_per_element_output.tsv.gz \\
            --per-guide-output perturbo_local_analysis_per_guide_output.tsv.gz \\
            --output-mudata inference_mudata.h5mu \\
            --v2-artifact-dir perturbo_v2_outputs \\
            --device ${params.INFERENCE_PERTURBO_DEVICE} \\
            --batch-size 0 \\
            --num-steps-control ${params.INFERENCE_PERTURBO_NUM_STEPS_CONTROL} \\
            --num-steps-betas ${params.INFERENCE_PERTURBO_NUM_STEPS_BETAS} \\
            --max-chunk-size ${params.INFERENCE_PERTURBO_LOCAL_MAX_CHUNK_CELLS} \\
            --perturbation-chunk-size 0 \\
            ${parallel_fits_arg} \\
            --element-gpu ${params.INFERENCE_PERTURBO_ELEMENT_GPU} \\
            --guide-gpu ${params.INFERENCE_PERTURBO_GUIDE_GPU} \\
            --jax-cache-dir ${params.INFERENCE_PERTURBO_JAX_CACHE_DIR} \\
            --size-factor-mode ${params.INFERENCE_PERTURBO_SIZE_FACTOR_MODE} \\
            --likelihood ${params.INFERENCE_PERTURBO_LIKELIHOOD} \\
            --prior ${params.INFERENCE_PERTURBO_PRIOR} \\
            ${save_model_params_arg}
        """
}
