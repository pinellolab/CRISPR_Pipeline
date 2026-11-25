process balance_mudata {
    cache 'lenient'
    cpus 2
    memory '128 GB'
    time '24h'
    publishDir "${params.outdir}/balanced_inputs", mode: 'copy'
    input:
    path source, stageAs: 'source/*'
    path parent, stageAs: 'parent/*'
    output:
    path 'sampled_full_mudata.h5mu', emit: full_mudata
    path 'balanced_inference_input.h5mu', emit: inference_input
    path 'downsampling_report.json', emit: report
    path 'selected_cells.tsv', emit: selected_cells
    script:
    """
    python ${projectDir}/bin/balance_mudata.py --source ${source} --parent-inference ${parent} \\
        --dataset '${params.BALANCED_dataset}' --cells ${params.BALANCED_cells} \\
        --rna-umis ${params.BALANCED_rna_umis} --seed ${params.BALANCED_seed}
    """
}
