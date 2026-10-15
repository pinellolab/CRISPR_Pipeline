nextflow.enable.dsl=2

include { seqSpecParser } from '../../../modules/local/seqSpecParser'
include { downloadReference } from '../../../modules/local/downloadReference'
include { mappingscRNA } from '../../../modules/local/mappingscRNA'
include { anndata_concat } from '../../../modules/local/anndata_concat'

workflow mapping_rna_pipeline {
    take:
    ch_rna
    ch_rna_seqspec
    ch_barcode_onlist
    parsed_covariate_file  // file: parsed covariate file from prepare_mapping_pipeline

    main:
    SeqSpecResult = seqSpecParser(
        ch_rna_seqspec,
        ch_barcode_onlist,
        'rna'
    )

    // A prebuilt index skips the download. The files are named as
    // downloadReference emits them, so the two branches are interchangeable.
    if (params.REFERENCE_index_local_path) {
        ref_dir = params.REFERENCE_index_local_path.replaceAll('/$', '')
        ref_transcriptome_idx = Channel.value(file("${ref_dir}/transcriptome_index.idx", checkIfExists: true))
        ref_t2g               = Channel.value(file("${ref_dir}/transcriptome_t2g.txt", checkIfExists: true))
        ref_cdna              = Channel.value(file("${ref_dir}/cdna.txt", checkIfExists: true))
        ref_nascent           = Channel.value(file("${ref_dir}/nascent_index.txt", checkIfExists: true))
    } else {
        DownloadRefResult = downloadReference(params.REFERENCE_transcriptome, params.use_igvf_reference)
        ref_transcriptome_idx = DownloadRefResult.transcriptome_idx
        ref_t2g               = DownloadRefResult.t2g_transcriptome_index
        ref_cdna              = DownloadRefResult.cdna
        ref_nascent           = DownloadRefResult.nascent_index
    }

    bc_replacement_ch = (params.replace_barcodes || params.rna_replace_barcodes) ? Channel.fromPath(params.bc_replacement_file).collect() : Channel.fromPath("${projectDir}/assets/dummy_bc_replacement.txt").collect()

    MappingOut = mappingscRNA(
        ch_rna,
        ref_transcriptome_idx,
        ref_t2g,
        ref_cdna,
        ref_nascent,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.barcode_file,
        bc_replacement_ch
    )

    ks_transcripts_out_dir_collected = MappingOut.ks_transcripts_out_dir
        .collect()
        .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }
    if (params.DEBUG_VAR) {
        ks_transcripts_out_dir_collected.view()
    }

    AnndataConcatenate = anndata_concat(
        parsed_covariate_file,
        ks_transcripts_out_dir_collected,
        params.rna_replace_barcodes
    )

    emit:
    trans_out_dir = MappingOut.ks_transcripts_out_dir
    ks_transcripts_out_dir_collected = MappingOut.ks_transcripts_out_dir
        .collect()
        .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }
    concat_anndata_rna = AnndataConcatenate.concat_anndata
    transcriptome_t2g = ref_t2g
}
