
process downloadGenome {
    input:
    val genome_path

    output:
    path "genome.fa.gz" , emit: genome

    script:
    """
        wget --continue --progress=dot:mega --tries=0 $genome_path -O genome.fa.gz
    """
}

