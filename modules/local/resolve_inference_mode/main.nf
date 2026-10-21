// The inference mode, settled once from the data for both methods.
//
// resolveControlGroupFromParams() (modules/local/control_group) validates the
// shared setting from the parameters alone and fails early on an unusable one.
// What it cannot see is the design: whether every cell carries one guide and
// whether any cell carries only controls. The PerTurbo adapter used to settle
// that on its own and override a declared MOI, while the SCEPTRE driver
// applied the same rule to the cis subset and never changed its control group,
// so the two could answer different questions on the same screen. This process
// measures the full guide matrix once and writes explicit values for both.
process resolve_inference_mode {
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
    path mudata            // every guide, the analysed cells: the MuData PerTurbo fits
    val setting            // INFERENCE_control_group, already validated
    val declared_moi       // Multiplicity_of_infection
    val perturbo_pool      // INFERENCE_PERTURBO_CRT_POOL (per-method override)
    val sceptre_group      // INFERENCE_SCEPTRE_control_group (per-method override)

    output:
    path "inference_mode.json", emit: json

    script:
    """
    resolve_inference_mode.py \\
        --mudata ${mudata} \\
        --setting '${setting}' \\
        --declared-moi '${declared_moi}' \\
        --perturbo-pool '${perturbo_pool}' \\
        --sceptre-group '${sceptre_group}' \\
        --output inference_mode.json
    """
}
