process FILTER_DEMO_SAMPLESHEET {
    tag "demo pre-run"
    container { params.containers.base }
    cpus 1
    memory '1 GB'

    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        return "${out.replaceAll('/$','')}/pipeline_info"
    }, mode: 'copy', overwrite: true
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        return out.replaceAll('/$','')
    }, mode: 'copy', overwrite: true, pattern: 'DEMO_MODE_WARNING.txt'

    input:
        path samplesheet
        val require_hash

    output:
        path "demo_samplesheet.*", emit: samplesheet
        path "DEMO_MODE_WARNING.txt", emit: warning

    script:
        def suffix = samplesheet.extension ?: 'csv'
        def hash_arg = require_hash ? '--require-hash' : ''
        """
        filter_demo_samplesheet.py \
            --input ${samplesheet} \
            --output "demo_samplesheet.${suffix}" \
            --warning-output DEMO_MODE_WARNING.txt \
            ${hash_arg}
        """
}
