
process skipGenomeDownload {
    cache 'lenient'
    debug true

    input:
    path genome_local_path
    
    output:
    path "genome.fa.gz", emit: genome

    script:
    """
    if [ -f "${genome_local_path}" ]; then
        if [ "${genome_local_path}" != "genome.fa.gz" ]; then
            cp "${genome_local_path}" genome.fa.gz
        else
            echo "File ${genome_local_path} found and renamed to genome.fa.gz"
        fi
    else
        echo "File ${genome_local_path} not found."
        exit 1
    fi
    """
}
