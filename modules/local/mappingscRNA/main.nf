process mappingscRNA {
    cache 'lenient'
    debug true
    stageOutMode 'copy'

    input:
    tuple val(meta), path(reads)
    path transcriptome_idx
    path transcriptome_t2g
    path cdna
    path nascent_idx
    path parsed_seqSpec_file
    path barcode_file
    path bc_replacement_file

    output:
    path "*_ks_transcripts_out", emit: ks_transcripts_out_dir
    path "*_ks_transcripts_out/counts_unfiltered*/adata.h5ad", emit: ks_transcripts_out_adata

    script:
    def batch = meta.measurement_sets
    def fastq_files = reads.join(' ')

    def workflow_args = params.scrna_workflow == "standard" ? 
        "--workflow standard" : 
        "--workflow nac -c1 ${cdna} -c2 ${nascent_idx}"

    def mm_flag = (params.use_multimapping) ? "--mm" : ""

    """
    echo "Processing ${batch} with ${fastq_files}"

    k_bin=\$(type -p kallisto)
    bustools_bin=\$(type -p bustools)
    chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})

    if [ -f "${bc_replacement_file}" ] && [ -s "${bc_replacement_file}" ]; then
        replacement_args="-r ${bc_replacement_file}"
    else
        replacement_args=""
    fi
    
    kb count \\
        -i ${transcriptome_idx} \\
        -g ${transcriptome_t2g} \\
        -o ${batch}_ks_transcripts_out \\
        -x \$chemistry \\
        -t ${task.cpus} \\
        ${workflow_args} \\
        ${mm_flag} \\
        --kallisto \$k_bin \\
        --bustools \$bustools_bin \\
        --overwrite \\
        --h5ad \\
        --sum total \\
        --verbose \\
        -w ${barcode_file} \\
        \$replacement_args \\
        ${fastq_files}

    echo "scRNA KB mapping Complete"
    """
}