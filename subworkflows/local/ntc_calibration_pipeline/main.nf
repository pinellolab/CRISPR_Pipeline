nextflow.enable.dsl=2

// The non-targeting calibration pass.
//
// PerTurbo already tests every non-targeting pseudo-element against every gene
// in its ordinary run. SCEPTRE tests only the pairs it is handed, so this
// subworkflow builds a second, small SCEPTRE input in which the non-targeting
// buckets are ordinary targets (under nt_cells, with the NT-only cells
// duplicated as the reference, mirroring PerTurbo's self-included control
// test), runs the same SCEPTRE processes on it, and joins the result to
// PerTurbo's global tables for the calibration and AUPRC report. The main
// run's pairs, families and tables are not touched.
//
// The SCEPTRE processes are included under their own names, not aliases, so
// every withName resource and container setting applies to this pass too.

include { prepare_ntc_calibration } from '../../../modules/local/prepare_ntc_calibration'
include { sceptre_chunk_prepare } from '../../../modules/local/sceptre_chunk_prepare'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { sceptre_chunk_merge } from '../../../modules/local/sceptre_chunk_merge'
include { ntc_calibration_report; ntc_calibration_skipped } from '../../../modules/local/ntc_calibration_report'

workflow ntc_calibration_pipeline {

    take:
    full_mudata                 // every guide: the MuData PerTurbo fits
    pairs_mudata                // carries uns['pairs_to_test'] (the main run's local pairs)
    mode                        // the resolved inference mode map (value channel)
    mode_json                   // inference_mode.json
    sceptre_local_per_guide     // the main run's SCEPTRE tables: AUPRC positives
    sceptre_local_per_element
    perturbo_global_per_guide   // PerTurbo's global tables: its non-targeting rows
    perturbo_global_per_element

    main:
    sort_paths = { paths -> paths.sort { a, b -> a.toString() <=> b.toString() } }

    Prep = prepare_ntc_calibration(full_mudata, pairs_mudata, mode)

    // The prepare step reports whether there is anything to calibrate; a screen
    // without controls skips the pass rather than failing the run.
    status = Prep.summary.map { f -> tuple(new groovy.json.JsonSlurper().parseText(f.text).status, f) }
    ok = status.filter { pair -> pair[0] == 'ok' }.map { pair -> pair[1] }
    skipped = status.filter { pair -> pair[0] != 'ok' }.map { pair -> pair[1] }

    // Every downstream input is gated through `ok`, so a skipped pass leaves
    // these channels empty and no SCEPTRE task is submitted.
    sceptre_input = Prep.sceptre_input.combine(ok).map { pair -> pair[0] }
    pairs_ok = Prep.pairs.combine(ok).map { pair -> pair[0] }

    Chunks = sceptre_chunk_prepare(sceptre_input)
    Results = inference_sceptre(Chunks.mudata_chunks.flatten(), mode)
    Merged = sceptre_chunk_merge(
        Results.per_guide_output.collect().map(sort_paths).combine(ok).map { pair -> pair[0] },
        Results.per_element_output.collect().map(sort_paths).combine(ok).map { pair -> pair[0] },
        sceptre_input,
        Chunks.chunk_manifest
    )

    Report = ntc_calibration_report(
        Merged.per_guide_output,
        Merged.per_element_output,
        perturbo_global_per_guide,
        perturbo_global_per_element,
        sceptre_local_per_guide,
        sceptre_local_per_element,
        pairs_ok,
        Prep.guide_targets,
        Prep.direct_targets,
        ok,
        mode_json,
        Prep.cell_counts_guide,
        Prep.cell_counts_element
    )
    Skipped = ntc_calibration_skipped(skipped, Prep.cell_counts_guide, Prep.cell_counts_element)

    emit:
    calibration_output = Report.calibration_output.mix(Skipped.calibration_output)
}
