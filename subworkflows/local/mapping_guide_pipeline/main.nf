nextflow.enable.dsl=2

include { seqSpecParser } from '../../../modules/local/seqSpecParser'
include { createGuideRef } from '../../../modules/local/createGuideRef'
include { mappingGuide } from '../../../modules/local/mappingGuide'
include { mappingGuideBaseEditing } from '../../../modules/local/mappingGuideBaseEditing'

include { anndata_concat } from '../../../modules/local/anndata_concat'

workflow mapping_guide_pipeline {
    take:
    ch_guide
    ch_guide_seqspec
    ch_barcode_onlist
    ch_guide_design
    parsed_covariate_file
    

    when:

    main:
    SeqSpecResult = seqSpecParser(
        ch_guide_seqspec,
        ch_barcode_onlist,
        'guide'
    )

    GuideRef = createGuideRef(ch_guide_design)

    if( params.is_BaseEditing ) {
            MappingOut = mappingGuideBaseEditing(
                ch_guide,
                SeqSpecResult.parsed_seqspec,
                SeqSpecResult.barcode_file,
                params.is_10x3v3,
                ch_guide_design

            )
        }
        else {
            MappingOut = mappingGuide(
                ch_guide,
                GuideRef.guide_index,
                GuideRef.t2g_guide,
                SeqSpecResult.parsed_seqspec,
                SeqSpecResult.barcode_file,
                params.is_10x3v3
            )
        }




    ks_guide_out_dir_collected = MappingOut.ks_guide_out_dir.collect()
    ks_guide_out_dir_collected.view()

    AnndataConcatenate = anndata_concat(
        parsed_covariate_file,
        ks_guide_out_dir_collected
    )

    emit:
    guide_out_dir = MappingOut.ks_guide_out_dir
    ks_guide_out_dir_collected = MappingOut.ks_guide_out_dir.collect()
    concat_anndata_guide = AnndataConcatenate.concat_anndata
}
