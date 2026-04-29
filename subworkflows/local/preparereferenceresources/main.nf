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
    def ch_reference_gtf

    if (params.transcriptome_index && params.transcriptome_t2g) {
        ch_index = channel.value(file(params.transcriptome_index, checkIfExists: true))
        ch_t2g = channel.value(file(params.transcriptome_t2g, checkIfExists: true))
        if (!params.reference_gtf) {
            error "`--reference_gtf` is required when using prebuilt transcriptome resources via `--transcriptome_index/--transcriptome_t2g`."
        }
        ch_reference_gtf = channel.value(file(params.reference_gtf, checkIfExists: true))
    } else {
        if (!params.reference_fasta || !params.reference_gtf) {
            error "Provide either `--transcriptome_index` and `--transcriptome_t2g` together, or provide `--reference_fasta` and `--reference_gtf` (or `--genome` that resolves them)."
        }
        KALLISTOBUSTOOLS_REF(
            ch_fasta,
            ch_gtf,
            ch_workflow_mode
        )
        ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_REF.out.versions.first())
        ch_index = KALLISTOBUSTOOLS_REF.out.index
        ch_t2g = KALLISTOBUSTOOLS_REF.out.t2g
        ch_reference_gtf = ch_gtf
    }

    emit:
    transcriptome_index = ch_index
    transcriptome_t2g   = ch_t2g
    reference_gtf       = ch_reference_gtf
    versions            = ch_versions
}
