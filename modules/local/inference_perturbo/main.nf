
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
    path "perturbo_cis_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_cis_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        def precomputed_element = params.INFERENCE_PERTURBO_PRECOMPUTED_PER_ELEMENT?.toString()?.trim()
        def precomputed_guide = params.INFERENCE_PERTURBO_PRECOMPUTED_PER_GUIDE?.toString()?.trim()
        def precomputed_mudata = params.INFERENCE_PERTURBO_PRECOMPUTED_MUDATA?.toString()?.trim()
        def supplied = [precomputed_element, precomputed_guide, precomputed_mudata].count { it }
        if (supplied != 0 && supplied != 3) {
            error "All three cis INFERENCE_PERTURBO_PRECOMPUTED_* paths must be set together"
        }
        def recovery_commands = precomputed_element ? """
        echo "Using precomputed cis PerTurbo results and MuData"
        cp --reflink=auto -f '${precomputed_element}' perturbo_cis_per_element_output.tsv.gz
        cp --reflink=auto -f '${precomputed_guide}' perturbo_cis_per_guide_output.tsv.gz
        cp --reflink=auto -f '${precomputed_mudata}' inference_mudata.h5mu
        """ : """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_cis_per_element_output.tsv.gz --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers 0 --efficiency_mode scaled --inference_type element

        # Run PerTurbo inference for per-guide results
        perturbo_inference.py ${mudata} perturbo_cis_per_guide_output.tsv.gz --batch_size ${params.INFERENCE_PERTURBO_BATCH_SIZE} --num_workers 0 --efficiency_mode scaled --inference_type guide

        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_results perturbo_cis_per_guide_output.tsv.gz \\
            --per_element_results perturbo_cis_per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
        """
        ${recovery_commands}
        """
}
