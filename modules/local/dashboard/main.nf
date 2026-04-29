process DASHBOARD {
    cache 'lenient'
    label 'process_single'

    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        return "${out}/pipeline_dashboard"
    }, mode: 'copy'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true, pattern: 'inference_mudata.h5mu'

    input:
    path guide_seqspec_plots
    path guide_fq_tbl
    path hashing_seqspec_plots
    path hashing_fq_tbl
    path mudata
    path gene_ann
    path gene_ann_filtered
    path guide_ann
    path hashing_ann
    path hashing_demux
    path hashing_unfiltered_demux
    path ks_transcripts_out_dir_collected
    path ks_guide_out_dir_collected
    path ks_hashing_out_dir_collected
    path additional_qc_dir
    path figures_dir
    path evaluation_output_dir
    path css
    path js
    path svg
    path controls_evaluation_output_dir
    path benchmark_output_dir

    output:
    tuple path("evaluation_output"), path("figures"), path("guide_seqSpec_plots"), path("hashing_seqSpec_plots"), path("additional_qc"), path("dashboard.html"), path("svg"), path("inference_mudata.h5mu"), path("benchmark_output")
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode = task.ext.mode ?: 'standard'
    def is_default = task.ext.default ?: false
    """
    set -euo pipefail

    # Normalize output folders that downstream expects.
    mkdir -p guide_seqSpec_plots hashing_seqSpec_plots figures evaluation_output svg additional_qc benchmark_output

    if [[ -d ${guide_seqspec_plots} ]]; then
        cp -rL ${guide_seqspec_plots}/* guide_seqSpec_plots/ 2>/dev/null || true
    fi
    if [[ -d ${hashing_seqspec_plots} ]]; then
        cp -rL ${hashing_seqspec_plots}/* hashing_seqSpec_plots/ 2>/dev/null || true
    fi
    if [[ -d ${figures_dir} ]]; then
        cp -rL ${figures_dir}/* figures/ 2>/dev/null || true
    fi
    if [[ -d ${evaluation_output_dir} ]]; then
        cp -rL ${evaluation_output_dir}/* evaluation_output/ 2>/dev/null || true
    fi
    if [[ -d ${controls_evaluation_output_dir} ]]; then
        cp -rL ${controls_evaluation_output_dir}/* evaluation_output/ 2>/dev/null || true
    fi
    if [[ -d ${svg} ]]; then
        cp -rL ${svg}/* svg/ 2>/dev/null || true
    fi
    if [[ -d ${additional_qc_dir} ]]; then
        cp -rL ${additional_qc_dir}/* additional_qc/ 2>/dev/null || true
    fi
    if [[ -d ${benchmark_output_dir} ]]; then
        cp -rL ${benchmark_output_dir}/* benchmark_output/ 2>/dev/null || true
    fi
    if [[ -f ${mudata} ]]; then
        cp -L ${mudata} inference_mudata.h5mu
    else
        touch inference_mudata.h5mu
    fi

    dashboard_entry.py \
        --mode ${mode} \
        --default ${is_default} \
        --output dashboard.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p guide_seqSpec_plots hashing_seqSpec_plots figures evaluation_output svg additional_qc benchmark_output
    touch dashboard.html
    touch inference_mudata.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
    END_VERSIONS
    """
}