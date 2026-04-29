include { KALLISTOBUSTOOLS_COUNT } from '../../../modules/nf-core/kallistobustools/count/main'

workflow QUANT_TRANSCRIPTOME {

    take:
    ch_reads        // channel: [ val(meta), [ reads ] ]  — one or more FASTQ files per sample
    ch_index        // channel: path  — kallisto index (.idx)
    ch_t2g          // channel: path  — transcript-to-gene mapping (t2g.txt)
    ch_whitelist    // channel: path or [] — barcode whitelist used by kb count
    ch_t1c          // channel: path  — cDNA t2c file (cdna_t2c.txt); pass [] to omit
    ch_t2c          // channel: path  — intron t2c file (intron_t2c.txt); pass [] to omit
    val_technology  // val: sequencing technology string, e.g. '10XV3'
    val_workflow    // val: kb workflow mode — 'standard', 'nac', or 'lamanno'
    val_prefix      // val: output directory prefix, e.g. 'A_ks_transcripts_out'

    main:
    def ch_versions = channel.empty()

    KALLISTOBUSTOOLS_COUNT (
        ch_reads,
        ch_index,
        ch_t2g,
        ch_whitelist,
        ch_t1c,
        ch_t2c,
        val_technology,
        val_workflow,
        val_prefix,
    )
    ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_COUNT.out.versions.first())

    emit:
    count    = KALLISTOBUSTOOLS_COUNT.out.count    // channel: [ val(meta), path("*.count") ]
    matrix   = KALLISTOBUSTOOLS_COUNT.out.matrix   // channel: path("*.count/*/*.mtx")
    versions = ch_versions                         // channel: path("versions.yml")
}
