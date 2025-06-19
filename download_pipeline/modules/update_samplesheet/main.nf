process update_samplesheet {
    cache 'lenient'
    debug true
    publishDir './pipeline_samplesheet'

    input:
    path per_sample_file

    output:
    path "updated_samplesheet.tsv", emit: updated_samplesheet

    script:
    """
    update_samplesheet.py --input ${per_sample_file} --output updated_samplesheet.tsv
    """
}
