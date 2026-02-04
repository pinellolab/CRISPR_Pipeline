process mappingGuide {
    tag "$meta.measurement_sets"
    cache 'lenient'
    debug true
    stageOutMode 'copy'

    input:
    tuple val(meta), path(reads)
    path guide_index
    path t2g_guide
    path parsed_seqSpec_file
    path barcode_file
    path bc_replacement_file
    val is_10xv3v   // should be "true" or "false"
    val spacer_tag

    output:
    path "*_ks_guide_out", emit: ks_guide_out_dir
    path "*_ks_guide_out/counts_unfiltered*/adata.h5ad", emit: ks_guide_out_adata

    script:
    def batch = meta.measurement_sets
    def fastq_files = reads.join(' ')

    // Check if spacer is valid (not null/empty and length > 1)
    def has_spacer  = (spacer_tag && spacer_tag.length() > 1) ? "true" : "false"
    """
    set -euo pipefail
    echo "Processing $batch"

    # 1. Determine the Base Chemistry
    # If it is 10xV3, we force that string. 
    # Otherwise, we extract it from the seqspec file.
    if [ "${is_10xv3v}" == "true" ]; then
        CHEM="10XV3"
        WORKFLOW="kite:10xFB"
    else
        RAW_CHEM=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})
        WORKFLOW="kite"
        
        # 2. Apply Spacer Logic (Modify Chemistry) if needed
        # Logic: If spacer exists, take the 3rd part of the string (transcript),
        # and set its 2nd and 3rd values to 0.
        # Ex: 0,0,16:0,16,26:1,20,30 -> 0,0,16:0,16,26:1,0,0
        
        if [ "${has_spacer}" == "true" ]; then
            echo "Spacer detected ('${spacer_tag}'). Modifying chemistry..."
            CHEM=\$(echo "\$RAW_CHEM" | awk -F: 'BEGIN{OFS=":"} {
                split(\$3, t, ","); 
                t[2]=0; 
                t[3]=0; 
                \$3=t[1]","t[2]","t[3]; 
                print 
            }')
        else
            CHEM="\$RAW_CHEM"
        fi
    fi

    echo "Final Chemistry: \$CHEM"

    # 3. Run kb count
    # We use the calculated \$CHEM and \$WORKFLOW variables

    if [ "${params.replace_barcodes}" = true ] && [ -f "${bc_replacement_file}" ] && [ -s "${bc_replacement_file}" ]; then
        replacement_args="-r ${bc_replacement_file}"
    else
        replacement_args=""
    fi

    kb count \\
        -i ${guide_index} \\
        -g ${t2g_guide} \\
        -o ${batch}_ks_guide_out \\
        -x \$CHEM \\
        -t ${task.cpus} \\
        --workflow \$WORKFLOW \\
        --h5ad \\
        --overwrite \\
        --verbose \\
        -w ${barcode_file} \\
        \$replacement_args \\
        ${fastq_files}

    echo "gRNA KB mapping Complete"
    """
}