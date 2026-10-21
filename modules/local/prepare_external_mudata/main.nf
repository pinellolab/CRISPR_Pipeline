process prepare_external_mudata {
    tag "${workflow.runName}"
    container { params.containers.base }
    publishDir { "${params.outdir}/external_mudata" }, mode: 'copy', overwrite: true

    input:
    path parent_manifest
    path parent_mudata, stageAs: 'parent/*'
    path edited_mudata, stageAs: 'edited/*'
    path gtf_reference
    path keep_cells, stageAs: 'filter/*'
    path covariates, stageAs: 'covariates/*'

    output:
    path 'analysis_input.h5mu', emit: mudata
    path 'external_mudata_report.json', emit: report

    script:
    def editedArg = params.ANALYSIS_input_mudata ? "--input-mudata ${edited_mudata}" : ''
    def cellsArg = params.ANALYSIS_keep_cells_tsv ? "--keep-cells ${keep_cells}" : ''
    def covariatesArg = params.ANALYSIS_covariates_tsv ? "--covariates ${covariates}" : ''
    """
    python ${projectDir}/bin/prepare_external_mudata.py \\
      --parent-manifest ${parent_manifest} --parent-mudata ${parent_mudata} \\
      --gtf ${gtf_reference} ${editedArg} ${cellsArg} ${covariatesArg}
    """
}
