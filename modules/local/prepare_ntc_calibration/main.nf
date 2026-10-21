// The SCEPTRE input for the non-targeting calibration pass: the NT buckets as
// ordinary targets, the three matched pair sets in uns['pairs_to_test'], and,
// under the nt_cells contrast, the NT-only cells duplicated once more as the
// reference population. The main run's pairs and tables are not touched.
process prepare_ntc_calibration {
    cache 'lenient'

    input:
        path full_mudata                         // every guide: the MuData PerTurbo fits
        path pairs_mudata, stageAs: 'pairs/*'    // carries uns['pairs_to_test'] (the local pairs)
        val mode                                 // the resolved inference mode map

    output:
        path "ntc_calibration_input.h5mu", optional: true, emit: sceptre_input
        path "ntc_calibration_pairs.tsv.gz", optional: true, emit: pairs
        path "ntc_guide_targets.tsv", emit: guide_targets
        path "ntc_direct_targets.tsv", emit: direct_targets
        path "ntc_cell_counts_per_guide.tsv", emit: cell_counts_guide
        path "ntc_cell_counts_per_element.tsv", emit: cell_counts_element
        path "ntc_calibration_design.json", emit: summary

    script:
        """
        prepare_ntc_calibration.py \\
            --full-mudata ${full_mudata} \\
            --pairs-mudata ${pairs_mudata} \\
            --sceptre-control-group ${mode.sceptre_control_group} \\
            --analysis-moi ${mode.analysis_moi} \\
            --max-local-pairs ${params.NTC_CALIBRATION_max_local_pairs} \\
            --max-global-pairs ${params.NTC_CALIBRATION_max_global_pairs} \\
            --seed ${params.NTC_CALIBRATION_seed} \\
            --output-mudata ntc_calibration_input.h5mu \\
            --pairs ntc_calibration_pairs.tsv.gz \\
            --guide-targets ntc_guide_targets.tsv \\
            --direct-targets ntc_direct_targets.tsv \\
            --cell-counts-per-guide ntc_cell_counts_per_guide.tsv \\
            --cell-counts-per-element ntc_cell_counts_per_element.tsv \\
            --summary ntc_calibration_design.json
        """
}
