process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2' :
        'quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  t2g
    path  whitelist
    path  t1c
    path  t2c
    val   technology
    val   workflow_mode
    val   prefix

    output:
    tuple val(meta), path ("${prefix}") , emit: count
    path "versions.yml"                 , emit: versions
    path "${prefix}/*/*.mtx"            , emit: matrix //Ensure that kallisto finished and produced outputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def cdna    = t1c ? "-c1 $t1c" : ''
    def intron  = t2c ? "-c2 $t2c" : ''
    def barcode = whitelist ? "-w $whitelist" : ''
    def memory  = task.memory.toGiga() - 1
    """
    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $barcode \\
        $cdna \\
        $intron \\
        -x $technology \\
        --workflow $workflow_mode \\
        $args \\
        -o ${prefix} \\
        -m ${memory}G \\
        ${reads.join( " " )}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}/counts_unfiltered/
    touch ${prefix}/counts_unfiltered/cells_x_genes.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
