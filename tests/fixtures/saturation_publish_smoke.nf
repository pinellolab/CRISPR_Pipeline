nextflow.enable.dsl = 2

include { sequencing_saturation } from '../../modules/local/sequencing_saturation'

workflow {
    sequencing_saturation(
        Channel.value(file("${projectDir}/assets/saturation_qc_empty")),
        Channel.value(file("${projectDir}/README.md")),
        Channel.value(file("${projectDir}/nextflow.config")),
        Channel.value(file("${projectDir}/CITATIONS.md"))
    )
}
