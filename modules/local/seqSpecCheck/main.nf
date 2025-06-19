process seqSpecCheck {
    cache 'lenient'
    debug true
    
    input:
    path(test_fastq_files_r1)
    path(test_fastq_files_r2)
    path(metadata)
    val data_type
    
    output:
    path "*_seqSpec_plots", emit: seqSpecCheck_plots
    path "*_position_table.csv", emit: position_table
    
    script:
    def r1_files = test_fastq_files_r1.join(' ')
    def r2_files = test_fastq_files_r2.join(' ')
    """
    echo "Checking fastq files for ${data_type}"
    seqSpecCheck.py --read1 ${r1_files} --read2 ${r2_files} --max_reads 100000 --metadata ${metadata} --plot
    mv seqSpec_plots ${data_type}_seqSpec_plots
    mv position_table.csv ${data_type}_position_table.csv
    """
}