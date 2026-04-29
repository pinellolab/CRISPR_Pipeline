// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process CREATEGUIDEREF {
    tag "$guide_metadata"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    path guide_metadata
    val  rev_comp
    val  spacer_tag

    output:
    path "guide_index.idx",   emit: guide_index
    path "t2guide.txt",       emit: t2g_guide
    path "guide_mismatch.fa", emit: guide_mismatch_fa
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ARGS=""
    PYTHON_SPACER=""
    if [ -n "${spacer_tag}" ]; then
        PYTHON_SPACER="--spacer ${spacer_tag}"
        ARGS="-k 31"
    fi

    guide_table.py --guide_table ${guide_metadata} --rev_comp ${rev_comp} \$PYTHON_SPACER

    kb ref \\
        -i guide_index.idx \\
        -f1 guide_mismatch.fa \\
        -g t2guide.txt \\
        \$ARGS \\
        --workflow kite \\
        guide_features.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch guide_index.idx
    touch t2guide.txt
    touch guide_mismatch.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
        pandas: "0"
        kallistobustools: "0"
    END_VERSIONS
    """
}
