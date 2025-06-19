
process downloadGTF {
    input:
    val gtf_url

    output:
    path "gencode_gtf.gtf.gz" , emit: gencode_gtf

    script:
    """
        wget --continue --progress=dot:mega --tries=0 $gtf_url -O gencode_gtf.gtf.gz
    """
}
