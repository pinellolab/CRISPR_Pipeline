
process inference_perturbo_trans {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, enabled: { params.INFERENCE_method == 'perturbo_trans' }

    input:
    path mudata
    val inference_method
    val dummy // fake dependency to force this to run after cis analysis

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_trans_per_element_output.parquet", emit: per_element_output
    path "perturbo_trans_per_guide_output.parquet", emit: per_guide_output
    path "perturbo_v2_outputs", optional: true, emit: perturbo_v2_outputs

    script:
        def save_model_params_arg = params.INFERENCE_PERTURBO_SAVE_MODEL_PARAMS ? '--save-model-params' : '--no-save-model-params'
        """
        perturbo_v2_pipeline_adapter.py \\
            --input ${mudata} \\
            --per-element-output perturbo_trans_per_element_output.parquet \\
            --per-guide-output perturbo_trans_per_guide_output.parquet \\
            --output-mudata inference_mudata.h5mu \\
            --v2-artifact-dir perturbo_v2_outputs \\
            --test-all-pairs \\
            --device ${params.INFERENCE_PERTURBO_DEVICE} \\
            --batch-size 0 \\
            --num-steps-control ${params.INFERENCE_PERTURBO_NUM_STEPS_CONTROL} \\
            --num-steps-betas ${params.INFERENCE_PERTURBO_NUM_STEPS_BETAS} \\
            --max-chunk-size ${params.INFERENCE_PERTURBO_MAX_CHUNK_CELLS} \\
            --perturbation-chunk-size ${params.INFERENCE_PERTURBO_PERTURBATION_CHUNK_SIZE} \\
            --size-factor-mode ${params.INFERENCE_PERTURBO_SIZE_FACTOR_MODE} \\
            --likelihood ${params.INFERENCE_PERTURBO_LIKELIHOOD} \\
            --prior ${params.INFERENCE_PERTURBO_PRIOR} \\
            ${save_model_params_arg}
        """
}
