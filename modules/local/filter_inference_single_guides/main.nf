process filter_inference_single_guides {
    cache 'lenient'
    cpus 2
    memory '64 GB'
    time '6h'
    container { params.containers.base }
    publishDir "${params.outdir}/inference_cell_filter", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/inference_cell_filter", mode: 'copy', pattern: '*.tsv'

    input:
    path mudata

    output:
    path 'single_guide_inference_input.h5mu', emit: mudata
    path 'inference_single_guide_filter.json', emit: report
    path 'inference_single_guide_cells.tsv', emit: cells

    script:
    """
    python ${projectDir}/bin/filter_inference_single_guides.py \\
        --input ${mudata} --dataset '${workflow.runName}' \\
        --min-umis ${params.INFERENCE_single_guide_min_umis} \\
        --min-ratio ${params.INFERENCE_single_guide_min_ratio}
    """
}
