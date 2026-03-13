
process inference_perturbo_trans {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }

    input:
    path mudata
    val inference_method
    val dummy // fake dependency to force this to run after cis analysis

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_trans_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_trans_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        PERTURBO_NUM_WORKERS=$(( ${task.cpus} > 1 ? ${task.cpus} - 1 : 0 ))

        # Run PerTurbo inference for per-element results
        perturbo_inference_chunked.py ${mudata} perturbo_trans_per_element_output.tsv.gz --chunk_size ${params.INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK} --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers ${PERTURBO_NUM_WORKERS} --efficiency_mode scaled --inference_type element --test_all_pairs
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference_chunked.py ${mudata} perturbo_trans_per_guide_output.tsv.gz --chunk_size ${params.INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK} --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers ${PERTURBO_NUM_WORKERS} --efficiency_mode scaled --inference_type guide --test_all_pairs
        
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_tsv perturbo_trans_per_guide_output.tsv.gz \\
            --per_element_tsv perturbo_trans_per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
}
