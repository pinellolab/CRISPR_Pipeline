process SEQSPECCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(reads)
    path  metadata         // seqspec YAML or guide-design TSV depending on data_type
    val   data_type        // string: 'rna', 'guide', or 'hash' — used as output prefix

    output:
    path "${data_type}_seqSpec_plots", emit: seqSpecCheck_plots
    path "${data_type}_position_table.csv", emit: position_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def max_reads  = task.ext.max_reads ?: 100000
    def read_pairs = reads.collate(2)
    def r1_files   = read_pairs.collect { it[0] }.join(' ')
    def r2_files   = read_pairs.findAll { it.size() > 1 }.collect { it[1] }.join(' ')
    """
    echo "Checking fastq files for ${data_type}"

    seqSpecCheck.py \\
        --read1 ${r1_files} \\
        --read2 ${r2_files} \\
        --max_reads ${max_reads} \\
        --metadata ${metadata} \\
        --plot \\
        ${args}

    mv seqSpec_plots ${data_type}_seqSpec_plots
    mv position_table.csv ${data_type}_position_table.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${data_type}_seqSpec_plots
    touch ${data_type}_seqSpec_plots/placeholder.png
    touch ${data_type}_position_table.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
