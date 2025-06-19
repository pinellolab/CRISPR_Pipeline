process download_portal {
    cache 'lenient'
    debug true

    input:
    path keypair_json
    val accession_id
    val download_option

    output:
    path "per-sample.tsv", emit: per_sample_file

    script:
    def download_args = ""
    def output_dir = ""

    if (download_option == "all") {
        download_args = ""
        output_dir = "--output-dir ./pipeline_inputs"
    } else if (download_option == "fastq") {
        download_args = "--file-types fastq"
        output_dir = "--output-dir ./fastq_dir"
    } else if (download_option == "other") {
        download_args = "--file-types other"
        output_dir = "--output-dir ./reference_files"
    } else {
        error "Invalid download_option: ${download_option}. Must be 'all', 'fastq', or 'other'"
    }

    """
    generate_per_sample.py --keypair ${keypair_json} --accession ${accession_id} --output per_sample.tsv
    python3 download_igvf.py --sample per_sample.tsv --keypair ${keypair_json} ${download_args} ${output_dir}
    """
}
