process preinference_checkpoint {
    tag "${workflow.runName}"
    container { params.containers.base }
    publishDir { "${params.outdir}/preinference_checkpoint/${workflow.runName}" }, mode: 'copy', overwrite: false

    input:
    path mudata_input
    path gtf_reference

    output:
    path 'preinference_mudata.h5mu', emit: mudata
    path 'preinference_manifest.json', emit: manifest

    script:
    """
    python ${projectDir}/bin/preinference_checkpoint.py \\
      --input ${mudata_input} --gtf ${gtf_reference} \\
      --run-name '${workflow.runName}' --commit '${workflow.commitId ?: workflow.manifest.version ?: 'unavailable'}'
    """
}
