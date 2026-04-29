process CREATEHASHINGREF {
    tag "$hashing_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path hashing_metadata

    output:
    path "hashing_index.idx",  emit: hashing_index
    path "t2g_hashing.txt",    emit: t2g_hashing
    path "hashing_mismatch.fa", emit: hashing_mismatch_fa
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    hashing_table.py --hashing_table ${hashing_metadata}

    kb ref \\
        -i hashing_index.idx \\
        -f1 hashing_mismatch.fa \\
        -g t2g_hashing.txt \\
        --workflow kite \\
        hashing_table.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch hashing_index.idx
    touch t2g_hashing.txt
    touch hashing_mismatch.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
        pandas: "0"
        kallistobustools: "0"
    END_VERSIONS
    """
}