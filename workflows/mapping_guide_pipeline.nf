nextflow.enable.dsl=2

include { seqSpecParser } from '../processes/seqSpecParser.nf'
include { createGuideRef } from '../processes/createGuideRef.nf'
include { mappingGuide } from '../processes/mappingGuide.nf'
include { anndata_concat } from '../processes/anndata_concat.nf'

workflow mapping_guide_pipeline{
    take:
    parsed_covariate_file
    genome

    main:
    SeqSpecResult = seqSpecParser(
        file("${params.seqspecs_directory}/${params.Guides_seqspec_yaml}"),
        file(params.seqspecs_directory),
        'guide')

    GuideRef = createGuideRef(file(params.guide_metadata))

    fastq_files_ch = Channel.fromPath(params.fastq_files_guide)
    batch_ch = Channel.fromList(params.batch)

    batch_fastq_ch = batch_ch.merge(fastq_files_ch)

    MappingOut = mappingGuide(
        batch_fastq_ch,
        GuideRef.guide_index,
        GuideRef.t2g_guide,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.barcode_file
        )
    
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
