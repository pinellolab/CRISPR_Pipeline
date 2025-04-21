
process createGuideRef {
    cache 'lenient'
    debug true

    input:
    path guide_metadata

    output:
    path "guide_index.idx" ,  emit: guide_index
    path "t2guide.txt" , emit: t2g_guide
    path "guide_mismatch.fa" , emit: guide_mismatch_fa

    script:

    """
        k_bin=\$(which kallisto)
        bustools_bin=\$(which bustools)
        guide_features_table=\$(guide_table.py --guide_table ${guide_metadata})
        kb ref -i guide_index.idx -f1 guide_mismatch.fa -g t2guide.txt --kallisto \$k_bin --bustools \$bustools_bin --workflow kite guide_features.txt
    """

}
