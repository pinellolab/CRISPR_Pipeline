nextflow.enable.dsl=2

include { seqSpecParser } from '../../../modules/local/seqSpecParser'
include { downloadGenome } from '../../../modules/local/downloadGenome'
include { createHashingRef } from '../../../modules/local/createHashingRef'
include { mappingHashing } from '../../../modules/local/mappingHashing'
include { anndata_concat } from '../../../modules/local/anndata_concat'

workflow mapping_hashing_pipeline {
    take:
    ch_hash
    ch_hash_seqspec
    ch_barcode_onlist
    ch_barcode_hashtag_map
    parsed_covariate_file

    main:
    SeqSpecResult = seqSpecParser(
        ch_hash_seqspec,
        ch_barcode_onlist,
        'hashing'
    )

    HashingRef = createHashingRef(ch_barcode_hashtag_map)

    MappingOut = mappingHashing(
        ch_hash,
        HashingRef.hashing_index,
        HashingRef.t2g_hashing,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.barcode_file
    )

    ks_hashing_out_dir_collected = MappingOut.ks_hashing_out_dir.collect()
    ks_hashing_out_dir_collected.view()

    AnndataConcatenate = anndata_concat(
        ks_hashing_out_dir_collected
        parsed_covariate_file,
        'hashing'
    )

    emit:
    hashing_out_dir = MappingOut.ks_hashing_out_dir
    ks_hashing_out_dir_collected = MappingOut.ks_hashing_out_dir.collect()
    concat_anndata_hashing = AnndataConcatenate.concat_anndata
}

