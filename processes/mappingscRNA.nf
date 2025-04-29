
process mappingscRNA {
    
    cache 'lenient'
    debug true

    input:
    each batch_fastq
    path transcriptome_idx
    path transcriptome_t2g
    path parsed_seqSpec_file
    path barcode_file

    output:
    path "*_ks_transcripts_out", emit: ks_transcripts_out_dir
    path "*_ks_transcripts_out/counts_unfiltered/adata.h5ad", emit: ks_transcripts_out_adata

    script:
        batch = batch_fastq[0]
        fastq_files = batch_fastq[1]
        """
        echo "Processing $batch with $fastq_files"

        k_bin=\$(type -p kallisto)
        bustools_bin=\$(type -p bustools)
        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})

        kb count -i ${transcriptome_idx} -g ${transcriptome_t2g} --verbose -w ${barcode_file} \\
                --h5ad --kallisto \$k_bin --bustools \$bustools_bin -x \$chemistry -o ${batch}_ks_transcripts_out -t ${task.cpus} \\
                ${fastq_files} --overwrite 

        echo "scRNA KB mapping Complete"
        """
}
