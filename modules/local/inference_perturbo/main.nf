
process inference_perturbo {
    cache 'lenient'
    publishDir './pipeline_outputs'

    input:
    path mudata
    val inference_method
    val efficiency_mode
    
    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv.gz", emit: per_element_output
    path "per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} per_element_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type element
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} per_guide_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type guide
        
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_tsv per_guide_output.tsv.gz \\
            --per_element_tsv per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
}
