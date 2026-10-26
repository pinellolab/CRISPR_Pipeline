
process mappingscRNA {

    cache 'lenient'
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
    path "*_ks_transcripts_out/counts_unfiltered${params.rna_replace_barcodes ? '_modified' : ''}/adata.h5ad", emit: ks_transcripts_out_adata

    script:
        def batch = meta.measurement_sets
        def fastq_files = reads.join(' ')
        def workflow_args = params.scrna_workflow == "standard" ? "--workflow standard" : "--workflow nac -c1 ${cdna} -c2 ${nascent_idx}"
        def mm_flag = (params.use_multimapping) ? "--mm" : ""
        // An empty scrna_strand leaves kb's own default in place.
        def strand_flag = params.scrna_strand ? "--strand ${params.scrna_strand}" : ""
        """
        echo "Processing $batch with $fastq_files"

        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})

        if [ "${params.rna_replace_barcodes}" = true ] && [ -f "${bc_replacement_file}" ] && [ -s "${bc_replacement_file}" ];
        then
                # bustools reads the replacement list with a plain ifstream, so
                # a gzipped file comes back as binary and stops with "Barcode
                # lengths not consistent in file".
                if gzip -t "${bc_replacement_file}" 2>/dev/null; then
                        gzip -dc "${bc_replacement_file}" > bc_replacement_list.txt
                        replacement_args="-r bc_replacement_list.txt"
                else
                        replacement_args="-r ${bc_replacement_file}"
                fi
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
                --overwrite \\
                --h5ad \\
                --sum total \\
                ${strand_flag} \\
                --verbose \\
                -w ${barcode_file} \\
                \$replacement_args \\
                ${fastq_files}

        echo "scRNA KB mapping Complete"
        """
}
