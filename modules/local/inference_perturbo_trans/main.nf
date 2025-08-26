
process inference_perturbo_trans {
    cache 'lenient'
    publishDir './pipeline_outputs'

    input:
    path mudata
    val inference_method
    val efficiency_mode

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv", emit: per_element_output
    path "per_guide_output.tsv", emit: per_guide_output


    script:
        """
        perturbo_inference.py ${mudata} per_element_output.tsv --test_all_pairs --mdata_output_fp inference_mudata.h5mu --efficiency_mode ${efficiency_mode} --inference_type element
        perturbo_inference.py ${mudata} per_guide_output.tsv --test_all_pairs --efficiency_mode ${efficiency_mode} --inference_type guide
        """
}
