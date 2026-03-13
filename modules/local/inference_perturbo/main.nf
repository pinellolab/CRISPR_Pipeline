
process inference_perturbo {
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
    
    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_results.tsv.gz", emit: perturbo_results
    path "perturbo_cis_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_cis_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        def perturboNumWorkers = task.cpus > 1 ? task.cpus - 1 : 0
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_cis_per_element_output.tsv.gz --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers ${perturboNumWorkers} --efficiency_mode scaled --inference_type element
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} perturbo_cis_per_guide_output.tsv.gz --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers ${perturboNumWorkers} --efficiency_mode scaled --inference_type guide
        
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_tsv perturbo_cis_per_guide_output.tsv.gz \\
            --per_element_tsv perturbo_cis_per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
}
