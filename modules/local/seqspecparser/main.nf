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

process SEQSPECPARSER {
    tag "$modality"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path seqspec_yaml
    path barcode_file
    val  modality

    output:
    path "${modality}_parsed_seqspec.txt", emit: parsed_seqspec
    path "barcode_file.txt",              emit: barcode_file
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parsing_guide_metadata.py \\
        --modalities ${modality} \\
        --yaml_file ${seqspec_yaml} \\
        --whitelist ${barcode_file} \\
        --output_file ${modality}_parsed_seqspec.txt

    cp ${barcode_file} barcode_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        seqspec: \$(seqspec --version 2>&1 | head -n 1 | sed 's/seqspec, version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${modality}_parsed_seqspec.txt
    cp ${barcode_file} barcode_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
        pandas: "0"
        seqspec: "0"
    END_VERSIONS
    """
}
