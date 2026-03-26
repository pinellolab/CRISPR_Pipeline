nextflow.enable.dsl=2

include { seqSpecCheck } from '../../../modules/local/seqSpecCheck'

workflow guideWorkflow {

    take:
    ch_guide
    ch_guide_design

    main:
    guide_seqSpecCheck = seqSpecCheck(
        ch_guide,
        ch_guide_design,
        'guide'
    )

    emit:
    guide_seqSpecCheck_plots = guide_seqSpecCheck.seqSpecCheck_plots
    guide_position_table = guide_seqSpecCheck.position_table
}

workflow seqSpecCheck_pipeline {
    take:
    ch_guide
    ch_guide_design

    main:
    guide = guideWorkflow(ch_guide, ch_guide_design)

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
}
