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
    val is_10xv3v
    val spacer_tag

    output:
    path "*_ks_guide_out"                              , emit: ks_guide_out_dir
    path "*_ks_guide_out/counts_unfiltered${params.replace_barcodes ? '_modified' : ''}/adata.h5ad" , emit: ks_guide_out_adata

    script:
    def batch = meta.measurement_sets
    def fastq_files = reads.join(' ')
    def counts_dir = "counts_unfiltered${params.replace_barcodes ? '_modified' : ''}"
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
        
        # 2. Apply Spacer Logic (Modify Chemistry) if needed.
        # If spacer exists, take the feature/guide part of the parsed chemistry,
        # preserve its read id, and set its start/end to 0.
        # Ex: 0,0,16:0,16,26:0,37,58 -> 0,0,16:0,16,26:0,0,0
        #     0,0,16:0,16,26:1,20,30 -> 0,0,16:0,16,26:1,0,0
        
        if [ "${has_spacer}" == "true" ]; then
            echo "Spacer tag detected ('${spacer_tag}'). Enabling automatic guide whole-read search."
            echo "The guide feature read is preserved from the parsed seqspec chemistry; only feature start/end are reset to 0,0."
            echo "Parsed guide chemistry: \$RAW_CHEM"
            CHEM=\$(echo "\$RAW_CHEM" | awk -F: 'BEGIN{OFS=":"} {
                if (NF >= 3) {
                    split(\$3, t, ",");
                    guide_read=t[1];
                    if (guide_read !~ /^[01]\$/) {
                        print "Could not infer guide feature read from parsed chemistry: " \$0 > "/dev/stderr";
                        exit 1;
                    }
                    \$3=guide_read",0,0";
                    print;
                    next;
                }
                if (NF == 2) {
                    n=split(\$2, t, ",");
                    if (n == 6) {
                        guide_read=t[4];
                        if (guide_read !~ /^[01]\$/) {
                            print "Could not infer guide feature read from parsed chemistry: " \$0 > "/dev/stderr";
                            exit 1;
                        }
                        \$2=t[1]","t[2]","t[3]":"guide_read",0,0";
                        print;
                        next;
                    }
                }
                print "Expected parsed chemistry in cb:umi:feature form, got: " \$0 > "/dev/stderr";
                exit 1;
            }')
            echo "Spacer-adjusted guide chemistry: \$CHEM"
        else
            CHEM="\$RAW_CHEM"
        fi
    fi

    echo "Final Chemistry: \$CHEM"

    if [ "${params.replace_barcodes}" = true ] && [ -f "${bc_replacement_file}" ] && [ -s "${bc_replacement_file}" ]; then
        replacement_args="-r ${bc_replacement_file}"
    else
        replacement_args=""
    fi

    # 3. Run kb count
    # We use the calculated \$CHEM and \$WORKFLOW variables
    kb count \\
        -i ${guide_index} \\
        -g ${t2g_guide} \\
        -w ${barcode_file} \\
        -o ${batch}_ks_guide_out \\
        -t ${task.cpus} \\
        -x \$CHEM \\
        --workflow \$WORKFLOW \\
        --h5ad \\
        --overwrite \\
        --verbose \\
        -w ${barcode_file} \\
        \${replacement_args} \\
        ${fastq_files}

    expected_h5ad="${batch}_ks_guide_out/${counts_dir}/adata.h5ad"
    if [ ! -s "\${expected_h5ad}" ]; then
        echo "ERROR: mappingGuide did not produce expected output: \${expected_h5ad}" >&2
        echo "Output tree for ${batch}_ks_guide_out:" >&2
        find "${batch}_ks_guide_out" -maxdepth 4 -type f -print >&2 || true
        exit 1
    fi
    """
}
