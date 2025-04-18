
process createHashingRef {
    cache 'lenient'

    input:
    path hashing_metadata

    output:
    path "hashing_index.idx" ,  emit: hashing_index
    path "t2g_hashing.txt" , emit: t2g_hashing
    path "hashing_mismatch.fa" , emit: hashing_mismatch_fa

    script:

    """
        k_bin=\$(which kallisto)
        bustools_bin=\$(which bustools)
        hashing_table=\$(hashing_table.py --hashing_table ${hashing_metadata})
        kb ref -i hashing_index.idx -f1 hashing_mismatch.fa -g t2g_hashing.txt --kallisto \$k_bin  --bustools \$bustools_bin --workflow kite hashing_table.txt
    """

}
