include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample'

workflow FASTQ_SUBSAMPLE {
    take:
    ch_fq_input

    main:

    FQ_SUBSAMPLE(ch_fq_input)

    ch_subsampled = FQ_SUBSAMPLE.out.fastq
        .map { meta, files ->
            def fastqs = files instanceof List ? files : [files]
            [meta, fastqs]
        }

    emit:
    subsample    = ch_subsampled
}
