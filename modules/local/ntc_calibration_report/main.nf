// The calibration and AUPRC report over matched non-targeting pairs, published
// as pipeline_outputs/calibration/. The SCEPTRE pass and the main run both name
// their tables sceptre_per_*_output.tsv.gz, so each pair is staged in its own
// directory.
process ntc_calibration_report {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true

    input:
        path sceptre_ntc_per_guide, stageAs: 'ntc_sceptre/*'
        path sceptre_ntc_per_element, stageAs: 'ntc_sceptre/*'
        path perturbo_global_per_guide
        path perturbo_global_per_element
        path sceptre_local_per_guide, stageAs: 'main_sceptre/*'
        path sceptre_local_per_element, stageAs: 'main_sceptre/*'
        path pairs
        path guide_targets
        path direct_targets
        path summary
        path inference_mode
        path cell_counts_guide
        path cell_counts_element

    output:
        path "calibration", emit: calibration_output

    script:
        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        export XDG_CACHE_HOME="\$PWD/.cache"
        mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME/fontconfig"
        ntc_calibration_report.py \\
            --sceptre-ntc-per-guide ${sceptre_ntc_per_guide} \\
            --sceptre-ntc-per-element ${sceptre_ntc_per_element} \\
            --perturbo-global-per-guide ${perturbo_global_per_guide} \\
            --perturbo-global-per-element ${perturbo_global_per_element} \\
            --sceptre-local-per-guide ${sceptre_local_per_guide} \\
            --sceptre-local-per-element ${sceptre_local_per_element} \\
            --pairs ${pairs} \\
            --guide-targets ${guide_targets} \\
            --direct-targets ${direct_targets} \\
            --summary-in ${summary} \\
            --inference-mode ${inference_mode} \\
            --cell-counts-per-guide ${cell_counts_guide} \\
            --cell-counts-per-element ${cell_counts_element} \\
            --outdir calibration
        """
}

// When the prepare step found nothing to calibrate (no non-targeting guides,
// or no control-only cell under the nt_cells contrast), the directory still
// exists, holding the reason and the cell-count sanity outputs.
process ntc_calibration_skipped {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true

    input:
        path summary
        path cell_counts_guide
        path cell_counts_element

    output:
        path "calibration", emit: calibration_output

    script:
        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        export XDG_CACHE_HOME="\$PWD/.cache"
        mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME/fontconfig"
        ntc_calibration_report.py --skipped \\
            --summary-in ${summary} \\
            --cell-counts-per-guide ${cell_counts_guide} \\
            --cell-counts-per-element ${cell_counts_element} \\
            --outdir calibration
        """
}
