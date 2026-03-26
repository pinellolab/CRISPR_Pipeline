process seqSpecCheck {
    cache 'lenient'
    debug true
    
    input:
    tuple val(meta), path(reads)
    path(metadata)
    val data_type
    
    output:
    path "*_seqSpec_plots", emit: seqSpecCheck_plots
    path "*_position_table.csv", emit: position_table
    
    script:
    def read_pairs = reads.collate(2)
    def r1_files = read_pairs.collect { it[0] }.join(' ')
    def r2_files = read_pairs.findAll { it.size() > 1 }.collect { it[1] }.join(' ')
    """
    echo "Checking fastq files for ${data_type}"
    seqSpecCheck.py --read1 ${r1_files} --read2 ${r2_files} --max_reads 100000 --metadata ${metadata} --plot
    mv seqSpec_plots ${data_type}_seqSpec_plots
    mv position_table.csv ${data_type}_position_table.csv
    """
}
