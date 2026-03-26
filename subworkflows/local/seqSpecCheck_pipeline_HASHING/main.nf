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

workflow hashWorkflow {
    take:
    ch_hash
    ch_barcode_hashtag_map

    main:
    hash_seqSpecCheck = seqSpecCheck(
        ch_hash,
        ch_barcode_hashtag_map,
        'hashing'
    )

    emit:
    hashing_seqSpecCheck_plots = hash_seqSpecCheck.seqSpecCheck_plots
    hashing_position_table = hash_seqSpecCheck.position_table
}


workflow seqSpecCheck_pipeline_HASHING {

    take:
    ch_guide
    ch_hash
    ch_guide_design
    ch_barcode_hashtag_map


    main:
    guide = guideWorkflow(ch_guide, ch_guide_design)
    hashing = hashWorkflow(ch_hash, ch_barcode_hashtag_map)

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
    hashing_seqSpecCheck_plots = hashing.hashing_seqSpecCheck_plots
    hashing_position_table = hashing.hashing_position_table
}
