//
// Subworkflow with functionality specific to the nf-core/crispr pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { FILTER_DEMO_SAMPLESHEET   } from '../../../modules/local/filter_demo_samplesheet'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    demo_mode         // boolean: Filter the samplesheet to one complete measurement set
    enable_hashing    // boolean: Require hash rows in the demo measurement set

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Copy the original samplesheet resolved from params/profile configuration
    //
    if (outdir && input) {
        copyOriginalSamplesheet(outdir, input)
    }

    ch_input_samplesheet = Channel.fromPath(input, checkIfExists: true)
    if (demo_mode) {
        log.warn """
================================================================================
DEMO MODE / PRE-RUN ONLY: only one complete measurement set will be processed.
Measurement set: complete set with the fewest scRNA FASTQ files
Do not report, publish, or treat these outputs as results from the full dataset.
================================================================================
""".stripIndent().trim()
        FILTER_DEMO_SAMPLESHEET(ch_input_samplesheet, enable_hashing)
        ch_samplesheet_file = FILTER_DEMO_SAMPLESHEET.out.samplesheet
    } else {
        ch_samplesheet_file = ch_input_samplesheet
    }

    //
    // Create channel from input file provided through params.input
    //

    samplesheet_sep = detectSamplesheetSeparator(input)

    ch_samplesheet_file
        .splitCsv(header: true, strip: true, sep: samplesheet_sep)
        .filter { row ->
            row.values().any { value -> value != null && value.toString().trim() }
        }
        .map { row ->
            validateSamplesheetRow(row)

            // Create a meta object with all required information
            def meta = [
                id: "sample_${row.file_modality}_${row.measurement_sets}",
                single_end: false,
                modality: row.file_modality,
                measurement_sets: row.measurement_sets,
                sequencing_run: row.sequencing_run,
                lane: row.lane,
                seqspec: row.seqspec,
                barcode_onlist: row.barcode_onlist,
                guide_design: row.guide_design,
                barcode_hashtag_map: row.barcode_hashtag_map
            ]

            def r1_path = row.R1_path ? file(row.R1_path) : null
            def r2_path = row.R2_path ? file(row.R2_path) : null

            // Return structured tuple with meta and file paths
            [ meta, r1_path, r2_path ]
        }
        .map {
            meta, r1_path, r2_path ->
                if (!r2_path) {
                    return [ meta.id, meta + [ single_end:true ], [ r1_path ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ r1_path, r2_path ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    demo_mode       // boolean: Mark completion as a reduced-data pre-run


    main:
    summary_params = [:]

    //
    // Completion email and summary
    //
    workflow.onComplete {
        def axiomSuccess = workflow?.success == true
        def axiomError = workflow?.errorMessage ?: 'unknown error'
        writeAxiomLifecycleEvent(
            'workflow_hook_complete',
            axiomSuccess ? 'SUCCEEDED' : 'FAILED',
            axiomSuccess ? 'Nextflow onComplete hook fired' : "Nextflow onComplete hook fired after failure: ${axiomError}"
        )
        if (demo_mode) {
            log.warn '''
================================================================================
DEMO MODE COMPLETED: these outputs were generated from one measurement set only.
They are pre-run diagnostics and must not be used as final full-dataset results.
See DEMO_MODE_WARNING.txt in the output directory.
================================================================================
'''.stripIndent().trim()
        }
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        writePipelineManifest(outdir)
        copyNextflowLog(outdir)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        def axiomError = workflow?.errorMessage ?: 'unknown error'
        writeAxiomLifecycleEvent(
            'workflow_hook_error',
            'FAILED',
            "Nextflow onError hook fired: ${axiomError}"
        )
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

// The sidecar tails this tiny file and performs network I/O outside Nextflow.
// Hook failures are swallowed so observability can never change pipeline status.
def writeAxiomLifecycleEvent(event_type, status, message) {
    if (!params.AXIOM_enabled || !params.AXIOM_event_file) {
        return null
    }
    try {
        def target = new File(params.AXIOM_event_file.toString())
        target.parentFile?.mkdirs()
        def payload = [
            _time: java.time.Instant.now().toString(),
            event_type: event_type,
            status: status,
            run_id: params.AXIOM_run_id ?: System.getenv('AXIOM_RUN_ID') ?: 'unknown',
            run_name: params.AXIOM_run_name ?: System.getenv('AXIOM_RUN_NAME') ?: workflow.runName,
            message: message.toString().take(2048)
        ]
        target << groovy.json.JsonOutput.toJson(payload) + System.lineSeparator()
    } catch (Exception error) {
        log.warn("Unable to write optional Axiom lifecycle event; pipeline continues: ${error.message}")
    }
}

//
// Copy the run log to the pipeline output directory
//
def copyNextflowLog(outdir) {
    def log_file = new File(workflow.launchDir.toString(), '.nextflow.log')
    if (!log_file.exists()) {
        log.warn("No .nextflow.log file found in ${workflow.launchDir}; skipping log export.")
        return null
    }

    nextflow.extension.FilesEx.copyTo(log_file.toPath(), "${outdir}/pipeline_info/nextflow.log")
}

//
// Write the exact pipeline source identity used by this completed run.
// This is a separate results artifact and never replaces the configuration
// supplied to Nextflow at launch.
//
def writePipelineManifest(outdir) {
    if (!outdir) {
        log.warn("No output directory configured; skipping pipeline manifest export.")
        return null
    }

    def project_dir = workflow.projectDir ?: workflow.launchDir
    def repository_path = project_dir.toString()

    def commit = getGitValue(repository_path, ['rev-parse', 'HEAD']) ?: workflow.commitId ?: workflow.manifest.version ?: 'unavailable'
    def branch = getGitValue(repository_path, ['symbolic-ref', '--short', '-q', 'HEAD']) ?: workflow.manifest.defaultBranch ?: 'unavailable'
    def repository_url = getGitValue(repository_path, ['remote', 'get-url', 'origin']) ?: workflow.manifest.homePage ?: 'unavailable'
    def version = commit == 'unavailable' ? 'unavailable' : commit.take(Math.min(12, commit.size()))
    def run_status = workflow.success ? 'success' : 'failed'
    def description = "Completed ${run_status} run ${workflow.runName}; Git commit ${commit}; project directory ${repository_path}"
    def content = """// Generated by CRISPR_Pipeline at workflow completion. Do not edit by hand.
// pipeline_command = '${escapeConfigString(workflow.commandLine ?: 'unavailable')}'
// nextflow_version = '${escapeConfigString(workflow.nextflow.version ?: 'unavailable')}'
// full_commit = '${escapeConfigString(commit)}'
// repository_path = '${escapeConfigString(repository_path)}'
// run_name = '${escapeConfigString(workflow.runName)}'
// completed_at = '${escapeConfigString(workflow.complete)}'
// success = ${workflow.success}
manifest {
    name = '${escapeConfigString(workflow.manifest.name ?: 'pinellolab/CRISPR_Pipeline')}'
    homePage = '${escapeConfigString(repository_url)}'
    defaultBranch = '${escapeConfigString(branch)}'
    version = '${escapeConfigString(version)}'
    description = '${escapeConfigString(description)}'
}
"""

    def target = new File("${outdir}/pipeline_info/pipeline_manifest.config")
    target.parentFile.mkdirs()
    target.text = content
    log.info("Wrote pipeline source manifest: ${target}")
}

def getGitValue(repository_path, arguments) {
    try {
        def command = ['git', '-C', repository_path.toString()] + arguments
        def process = new ProcessBuilder(command).redirectErrorStream(true).start()
        if (!process.waitFor(20, java.util.concurrent.TimeUnit.SECONDS) || process.exitValue() != 0) {
            return null
        }
        def value = process.inputStream.getText('UTF-8').trim()
        return value ? value : null
    } catch (Exception ignored) {
        return null
    }
}

def escapeConfigString(value) {
    return value.toString().replace('\\', '\\\\').replace("'", "\\'").replaceAll(/[\r\n]+/, ' ')
}

//
// Copy the original samplesheet after profile/config resolution
//
def copyOriginalSamplesheet(outdir, input) {
    def samplesheet = file(input)
    if (!samplesheet.exists()) {
        log.warn("Original samplesheet '${input}' was not found; skipping samplesheet export.")
        return null
    }

    def samplesheet_path = samplesheet instanceof java.nio.file.Path ? samplesheet : samplesheet.toPath()
    def name = samplesheet_path.getFileName().toString()
    def dot_index = name.lastIndexOf('.')
    def suffix = dot_index >= 0 ? name.substring(dot_index) : ''
    nextflow.extension.FilesEx.copyTo(samplesheet_path, "${outdir}/pipeline_info/original_samplesheet${suffix}")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size() == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Path-list order contributes to the Nextflow task hash and rendered
    // mapping command. Sort whole read pairs (never individual R1/R2 paths)
    // so a reordered/regenerated samplesheet still reuses the mapping cache.
    def entries = (0..<metas.size()).collect { index ->
        [meta: metas[index], reads: fastqs[index]]
    }.sort { left, right ->
        def left_key = [
            left.meta.sequencing_run?.toString() ?: '',
            left.meta.lane?.toString() ?: '',
            left.reads.collect { it.toString() }.join('\t')
        ].join('\t')
        def right_key = [
            right.meta.sequencing_run?.toString() ?: '',
            right.meta.lane?.toString() ?: '',
            right.reads.collect { it.toString() }.join('\t')
        ].join('\t')
        left_key <=> right_key
    }

    return [ entries[0].meta, entries.collect { it.reads } ]
}

//
// Detect whether the samplesheet is comma- or tab-separated.
//
def detectSamplesheetSeparator(input) {
    def input_file = file(input)
    def header = input_file.text.readLines().find { line -> line?.trim() }

    if (header?.contains('\t')) {
        return '\t'
    }
    return ','
}

//
// Validate required samplesheet columns before constructing metadata.
//
def validateSamplesheetRow(row) {
    def required_columns = [
        'R1_path',
        'file_modality',
        'measurement_sets',
        'sequencing_run',
        'lane',
        'seqspec',
        'barcode_onlist',
        'guide_design'
    ]
    // lane is retained in the input schema but may be empty; samples are grouped by modality and measurement set.
    def required_value_columns = required_columns - ['lane']
    def missing_columns = required_columns.findAll { column -> !row.containsKey(column) }
    if (missing_columns) {
        error("Input samplesheet is missing required column(s): ${missing_columns.join(', ')}. Check that the file is comma- or tab-delimited and has the expected header.")
    }

    def missing_values = required_value_columns.findAll { column -> !row[column]?.toString()?.trim() }
    if (missing_values) {
        def row_context = [
            "file_modality=${row.file_modality ?: '<empty>'}",
            "measurement_sets=${row.measurement_sets ?: '<empty>'}",
            "sequencing_run=${row.sequencing_run ?: '<empty>'}"
        ].join(', ')
        error("Input samplesheet row is missing required value(s): ${missing_values.join(', ')}. ${row_context}")
    }

    return true
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",


            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [


        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
