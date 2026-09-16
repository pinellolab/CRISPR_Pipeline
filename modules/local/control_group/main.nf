// One setting for the cells a perturbation is compared against.
//
// SCEPTRE and PerTurbo make the same statistical choice under different names.
// params.INFERENCE_control_group states it once, in SCEPTRE's vocabulary, and
// this is the mapping the pipeline hands to the two processes:
//
//     setting        SCEPTRE control_group   PerTurbo --crt-pool
//     -----------    ---------------------   -------------------
//     nt_cells       nt_cells                control-anchored
//     complement     complement              all-cells
//     auto           from the declared MOI   from the declared MOI
//
// 'auto' reads params.Multiplicity_of_infection: 'low' means the non-targeting
// cells are the reference population, 'high' means the complement of the
// perturbation is (the only contrast SCEPTRE offers at high MOI). When the MOI
// is neither, PerTurbo measures the design from the data and SCEPTRE falls back
// to 'complement', which is valid for every design.
//
// This mirrors bin/control_group.py, which is the version the PerTurbo adapter
// imports and the version the unit tests exercise; the two tables are checked
// against each other in tests/test_control_group_resolution.py. Change both.

def controlGroupSettings() {
    return ['auto', 'nt_cells', 'complement']
}

// The historical per-method defaults. A per-method parameter still sitting on
// its historical value is "not set" for precedence: every config written before
// the shared setting existed carries these.
def controlGroupPerturboHistoricalDefault() { return 'from-moi' }
def controlGroupSceptreHistoricalDefault() { return 'complement' }

def controlGroupSceptreBySetting() { return ['nt_cells': 'nt_cells', 'complement': 'complement'] }
def controlGroupPerturboBySetting() { return ['nt_cells': 'control-anchored', 'complement': 'all-cells'] }
def controlGroupPerturboByMoi() { return ['low': 'control-anchored', 'high': 'all-cells'] }
def controlGroupSceptreByMoi() { return ['low': 'nt_cells', 'high': 'complement'] }

def normalizeMoi(moi) {
    def text = moi == null ? '' : moi.toString().trim().toLowerCase()
    return (text in ['low', 'high']) ? text : 'unknown'
}

def isControlGroupOverride(value, historicalDefault) {
    if (value == null) {
        return false
    }
    def text = value.toString().trim()
    return text && text.toLowerCase() != historicalDefault.toLowerCase()
}

// Resolve the shared setting into one value per method. Returns a Map with
// sceptre_control_group, perturbo_crt_pool, the provenance of each, and the one
// log line that says what happened and why.
def resolveControlGroup(setting, moi, perturboPool = null, sceptreGroup = null) {
    def resolvedSetting = setting == null ? 'auto' : setting.toString().trim().toLowerCase()
    if (!(resolvedSetting in controlGroupSettings())) {
        error(
            "INFERENCE_control_group='${setting}' is not a control group. " +
            "Accepted values: ${controlGroupSettings().join(', ')}. " +
            "'nt_cells' compares each perturbation with the non-targeting cells, " +
            "'complement' with every other cell, and 'auto' takes whichever the " +
            "declared Multiplicity_of_infection implies."
        )
    }
    def resolvedMoi = normalizeMoi(moi)

    def sceptre
    def perturbo
    def perturboEffective
    def reason
    if (resolvedSetting == 'auto') {
        sceptre = controlGroupSceptreByMoi().getOrDefault(resolvedMoi, 'complement')
        // 'from-moi' hands the mapping to the adapter, which also gets to override
        // it from the assignments it can see. Naming a pool here would switch that off.
        perturbo = controlGroupPerturboHistoricalDefault()
        perturboEffective = controlGroupPerturboByMoi().getOrDefault(resolvedMoi, 'auto')
        reason = resolvedMoi == 'unknown'
            ? "auto with no declared low/high MOI: PerTurbo measures the design, SCEPTRE falls back to 'complement', the only contrast valid for every design"
            : "auto from declared MOI '${resolvedMoi}'"
    }
    else {
        sceptre = controlGroupSceptreBySetting()[resolvedSetting]
        perturbo = controlGroupPerturboBySetting()[resolvedSetting]
        perturboEffective = perturbo
        reason = "INFERENCE_control_group='${resolvedSetting}' set explicitly"
    }

    def sceptreProvenance = resolvedSetting == 'auto' ? 'auto' : 'explicit'
    def perturboProvenance = sceptreProvenance
    def overrides = [:]
    if (isControlGroupOverride(sceptreGroup, controlGroupSceptreHistoricalDefault())) {
        sceptre = sceptreGroup.toString().trim()
        sceptreProvenance = 'per-method-override'
        overrides['INFERENCE_SCEPTRE_control_group'] = sceptre
    }
    if (isControlGroupOverride(perturboPool, controlGroupPerturboHistoricalDefault())) {
        perturbo = perturboPool.toString().trim()
        perturboEffective = perturbo
        perturboProvenance = 'per-method-override'
        overrides['INFERENCE_PERTURBO_CRT_POOL'] = perturbo
    }

    // The one combination no amount of configuration can deliver: SCEPTRE has no
    // non-targeting-cell contrast at high MOI. Fail here rather than let the R
    // driver quietly substitute the complement, which is how the two methods
    // ended up answering different questions.
    if (sceptre == 'nt_cells' && resolvedMoi == 'high') {
        error(
            "control group 'nt_cells' is not available for a high-MOI screen: SCEPTRE's " +
            "non-targeting-cell contrast needs one perturbation per cell. Declared " +
            "Multiplicity_of_infection='high'. Use INFERENCE_control_group='complement', or " +
            "declare the screen low-MOI if that is what it is."
        )
    }

    def perturboDisplay = perturbo == perturboEffective ? perturbo : "${perturbo} -> ${perturboEffective}"
    def line = "Control group: declared MOI '${resolvedMoi}', INFERENCE_control_group='${resolvedSetting}'" +
        " -> SCEPTRE control_group='${sceptre}', PerTurbo --crt-pool='${perturboDisplay}'" +
        " (${reason})."
    if (overrides) {
        def which = overrides.sort { entry -> entry.key }.collect { k, v -> "${k}='${v}'" }.join(', ')
        line += " ${which} overrides the shared setting for that method only:" +
            " SCEPTRE and PerTurbo are now deliberately inconsistent and their calls" +
            " are not comparable."
    }
    else {
        line += ' Both methods contrast against the same cells.'
    }

    return [
        setting: resolvedSetting,
        moi: resolvedMoi,
        sceptre_control_group: sceptre,
        perturbo_crt_pool: perturbo,
        perturbo_pool_effective: perturboEffective,
        sceptre_provenance: sceptreProvenance,
        perturbo_provenance: perturboProvenance,
        overrides: overrides,
        reason: reason,
        log_line: line,
    ]
}

// The pipeline's own resolution, from params. Call it once per run.
def resolveControlGroupFromParams() {
    return resolveControlGroup(
        params.INFERENCE_control_group,
        params.Multiplicity_of_infection,
        params.INFERENCE_PERTURBO_CRT_POOL,
        params.INFERENCE_SCEPTRE_control_group
    )
}
