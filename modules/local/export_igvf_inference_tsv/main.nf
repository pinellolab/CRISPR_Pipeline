process export_igvf_inference_tsv {
    tag 'four canonical inference tables'
    publishDir "${params.outdir}/results/pipeline_results", mode: 'move', overwrite: true

    input:
    path local_element
    path local_guide
    path global_element
    path global_guide

    output:
    path '*_output.tsv', emit: igvf_tables

    script:
    """
    export_igvf_inference_tsv.py --sources ${local_element} ${local_guide} ${global_element} ${global_guide} --output-dir .
    """
}
