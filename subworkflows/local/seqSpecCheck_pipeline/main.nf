nextflow.enable.dsl=2

include { seqSpecCheck } from '../../../modules/local/seqSpecCheck'

workflow guideWorkflow {

    take:
    ch_guide

    main:
    ch_guide_split = ch_guide.map { meta, reads ->
        def r1 = reads[0]  // First file is R1
        def r2 = reads[1]  // Second file is R2
        [meta, r1, r2]
    }

    guide_seqSpecCheck = seqSpecCheck(
        ch_guide_split.map { meta, r1, r2 -> r1 },  // R1 channel
        ch_guide_split.map { meta, r1, r2 -> r2 },  // R2 channel
        file(params.METADATA_sgRNA),
        'guide'
    )

    emit:
    guide_seqSpecCheck_plots = guide_seqSpecCheck.seqSpecCheck_plots
    guide_position_table = guide_seqSpecCheck.position_table
}

workflow seqSpecCheck_pipeline {
    take:
    ch_guide

    main:
    guide = guideWorkflow(ch_guide)

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
}
