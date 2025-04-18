nextflow.enable.dsl=2

include { seqSpecParser } from '../processes/seqSpecParser.nf'
include { downloadGenome } from '../processes/downloadGenome.nf'
include { createHashingRef } from '../processes/createHashingRef.nf'
include { mappingHashing } from '../processes/mappingHashing.nf'
include { anndata_concat } from '../processes/anndata_concat.nf'

workflow mapping_hashing_pipeline{
    take:
    parsed_covariate_file
    genome

    main:
    SeqSpecResult = seqSpecParser(
        file("${params.seqspecs_directory}/${params.Hash_seqspec_yaml}"),
        file(params.seqspecs_directory),
        'hashing')

    HashingRef = createHashingRef(file(params.hashing_metadata))

    fastq_files_ch = Channel.fromPath(params.fastq_files_hashing)
    batch_ch = Channel.fromList(params.batch)

    batch_fastq_ch = batch_ch.merge(fastq_files_ch)

    MappingOut = mappingHashing(
        batch_fastq_ch,
        HashingRef.hashing_index,
        HashingRef.t2g_hashing,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.barcode_file
        )
    
    ks_hashing_out_dir_collected = MappingOut.ks_hashing_out_dir.collect()
    ks_hashing_out_dir_collected.view()
    
    AnndataConcatenate = anndata_concat(
        parsed_covariate_file,
        ks_hashing_out_dir_collected
    )

    emit:
    hashing_out_dir = MappingOut.ks_hashing_out_dir
    ks_hashing_out_dir_collected = MappingOut.ks_hashing_out_dir.collect()
    concat_anndata_hashing = AnndataConcatenate.concat_anndata
}
