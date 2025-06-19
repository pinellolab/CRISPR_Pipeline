
process skipGTFDownload {
    cache 'lenient'
    debug true

    input:
    path gtf_local_path
    
    output:
    path "gencode_gtf.gtf.gz", emit: gtf

    script:
    """
    if [ -f "${gtf_local_path}" ]; then
        if [ "${gtf_local_path}" != "gencode_gtf.gtf.gz" ]; then
            cp "${gtf_local_path}" gencode_gtf.gtf.gz
        else
            echo "File ${gtf_local_path} found and renamed to gencode_gtf.gtf.gz"
        fi
    else
        echo "File ${gtf_local_path} not found."
        exit 1
    fi
    """
}
