include { KALLISTOBUSTOOLS_REF } from '../../../modules/nf-core/kallistobustools/ref/main'

workflow PREPAREREFERENCERESOURCES {

    take:
    ch_fasta
    ch_gtf
    ch_workflow_mode

    main:
    def ch_versions = channel.empty()
    def ch_index
    def ch_t2g

    if (params.transcriptome_index && params.transcriptome_t2g) {
        ch_index = channel.value(file(params.transcriptome_index, checkIfExists: true))
        ch_t2g = channel.value(file(params.transcriptome_t2g, checkIfExists: true))
    } else {
        KALLISTOBUSTOOLS_REF(
            ch_fasta,
            ch_gtf,
            ch_workflow_mode
        )
        ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_REF.out.versions.first())
        ch_index = KALLISTOBUSTOOLS_REF.out.index
        ch_t2g = KALLISTOBUSTOOLS_REF.out.t2g
    }

    emit:
    transcriptome_index = ch_index
    transcriptome_t2g   = ch_t2g
    reference_gtf       = ch_gtf
    versions            = ch_versions
}
