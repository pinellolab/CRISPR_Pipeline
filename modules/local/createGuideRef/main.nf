process createGuideRef {
    debug true

    input:
    path guide_metadata
    val rev_comp
    val spacer_tag

    output:
    path "guide_index.idx"   , emit: guide_index
    path "t2guide.txt"       , emit: t2g_guide
    path "guide_mismatch.fa" , emit: guide_mismatch_fa
    
    script:
    """
    # Define optional flag based on spacer_tag availability
    ARGS=""
    PYTHON_SPACER=""
    echo "Spacer tag provided: '${spacer_tag}'"
    if [ -z "${spacer_tag}" ]; then
        echo ' No spacer tag and keeping kmer size automatic' 
    else
        PYTHON_SPACER="--spacer ${spacer_tag}"
        ARGS="-k 31"
    fi
    
    echo "Spacer tag being used: '${spacer_tag}'"

    # Run python script
    guide_table.py  --guide_table ${guide_metadata} --rev_comp ${rev_comp} \$PYTHON_SPACER 
    echo 'checking head of guide features file'
    head guide_features.txt
    echo 'args used'
    echo \$ARGS
    echo "Guide features table created."


    # Run kb ref
    kb ref \\
        -i guide_index.idx \\
        -f1 guide_mismatch.fa \\
        -g t2guide.txt \\
        \$ARGS \\
        --workflow kite \\
        guide_features.txt
    """
}