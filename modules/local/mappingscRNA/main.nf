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

    output:
    path "*_ks_transcripts_out", emit: ks_transcripts_out_dir
    path "*_ks_transcripts_out/counts_unfiltered/adata.h5ad", emit: ks_transcripts_out_adata

    script:
    def batch = meta.measurement_sets
    def fastq_files = reads.join(' ')

    def workflow_args = params.scrna_workflow == "standard" ? 
        "--workflow standard" : 
        "--workflow nac -c1 ${cdna} -c2 ${nascent_idx}"
    
    def replacement_args = (params.replace_barcodes && params.bc_replacement_file != '') ?
        "-r ${params.bc_replacement_file}" : ""

    def mm_flag = (params.scrna_workflow == "standard" && params.use_multimapping) ? "--mm" : ""
    // If nascent workflow, always use --mm
    if (params.scrna_workflow == "nac") {
        mm_flag = "--mm"
    }
    
    """
    echo "Processing ${batch} with ${fastq_files}"
    k_bin=\$(type -p kallisto)
    bustools_bin=\$(type -p bustools)
    chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})
    
    kb count \\
        -i ${transcriptome_idx} \\
        -g ${transcriptome_t2g} \\
        ${workflow_args} \\
        ${mm_flag} \\
        ${replacement_args} \\
        -x \$chemistry \\
        -w ${barcode_file} \\
        -o ${batch}_ks_transcripts_out \\
        -t ${task.cpus} \\
        ${fastq_files} \\
        --h5ad \\
        --kallisto \$k_bin \\
        --bustools \$bustools_bin \\
        --overwrite \\
        --verbose

    echo "scRNA KB mapping Complete"
    """
}