nextflow.enable.dsl=2

include { download_portal } from './modules/download_portal'
include { download_portal_no_seqspec } from './modules/download_portal_no_seqspec'
include { update_samplesheet } from './modules/update_samplesheet'

workflow {

    main:
    if( params.SEQSPEC_ON_PORTAL == "true" ){
        DownloadPortal = download_portal(
            file(params.keypair_json),
            params.accession_id,
            params.download_option
            )
    } else {
        DownloadPortal = download_portal_no_seqspec(
            file(params.keypair_json),
            params.accession_id,
            params.download_option,
            file(params.SEQUENCE_PARSING_hash_seqspec_yaml),
            file(params.SEQUENCE_PARSING_scRNA_seqspec_yaml),
            file(params.SEQUENCE_PARSING_sgRNA_seqspec_yaml)
            )
    }

    update_samplesheet(DownloadPortal.per_sample_file)
}
