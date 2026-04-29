include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample'

workflow FASTQ_SUBSAMPLE {
    take:
    ch_inputs

    main:

    ch_fq_input = ch_inputs
        .map { meta, r1_path, r2_path ->
            def reads = [r1_path, r2_path]
            [meta, reads]
        }

    FQ_SUBSAMPLE(ch_fq_input)

    ch_subsampled = FQ_SUBSAMPLE.out.fastq
        .map { meta, files ->
            def r1_path = files instanceof List ? files[0] : files
            def r2_path = files instanceof List && files.size() > 1 ? files[1] : null
            [meta, r1_path, r2_path]
        }

    emit:
    reads    = ch_subsampled
}
