
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
    path "perturbo_trans_per_element_output.*", emit: per_element_output
    path "perturbo_trans_per_guide_output.*", emit: per_guide_output

    script:
        def results_ext = params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT == 'parquet' ? 'parquet' : 'tsv.gz'
        def precomputed_element = params.INFERENCE_PERTURBO_TRANS_PRECOMPUTED_PER_ELEMENT?.toString()?.trim()
        def precomputed_guide = params.INFERENCE_PERTURBO_TRANS_PRECOMPUTED_PER_GUIDE?.toString()?.trim()
        def precomputed_mudata = params.INFERENCE_PERTURBO_TRANS_PRECOMPUTED_MUDATA?.toString()?.trim()
        if ((precomputed_element && !precomputed_guide) || (!precomputed_element && precomputed_guide)) {
            error "Both INFERENCE_PERTURBO_TRANS_PRECOMPUTED_PER_ELEMENT and INFERENCE_PERTURBO_TRANS_PRECOMPUTED_PER_GUIDE must be set together"
        }
        if (precomputed_mudata && !precomputed_element) {
            error "INFERENCE_PERTURBO_TRANS_PRECOMPUTED_MUDATA requires both precomputed trans result tables"
        }
        if (precomputed_element) {
            def expected_suffix = ".${results_ext}"
            if (!precomputed_element.endsWith(expected_suffix) || !precomputed_guide.endsWith(expected_suffix)) {
                error "Precomputed trans PerTurbo result suffixes must match INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT='${params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT}'"
            }
        }
        def result_commands = precomputed_element ? """
        echo "Using precomputed trans PerTurbo results"
        cp -f '${precomputed_element}' perturbo_trans_per_element_output.${results_ext}
        cp -f '${precomputed_guide}' perturbo_trans_per_guide_output.${results_ext}
        """ : """
        # Run PerTurbo inference for per-element results
        perturbo_inference_chunked.py ${mudata} perturbo_trans_per_element_output.${results_ext} --chunk_size ${params.INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK} --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers 0 --efficiency_mode scaled --inference_type element --test_all_pairs
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference_chunked.py ${mudata} perturbo_trans_per_guide_output.${results_ext} --chunk_size ${params.INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK} --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers 0 --efficiency_mode scaled --inference_type guide --test_all_pairs
        """
        def mudata_command = precomputed_mudata ? """
        echo "Using precomputed trans PerTurbo MuData"
        cp --reflink=auto -f '${precomputed_mudata}' inference_mudata.h5mu
        """ : """
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_results perturbo_trans_per_guide_output.${results_ext} \\
            --per_element_results perturbo_trans_per_element_output.${results_ext} \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
        """
        ${result_commands}
        ${mudata_command}
        """
}
