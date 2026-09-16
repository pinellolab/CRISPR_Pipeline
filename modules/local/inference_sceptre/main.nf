
process inference_sceptre {

    input:
    path mudata_fp
    // The cells a perturbation is compared against: 'nt_cells' or 'complement',
    // resolved once from params.INFERENCE_control_group in the inference
    // subworkflow (modules/local/control_group) so PerTurbo's CRT pool and this
    // are the same choice. The R driver honours it and refuses the one
    // combination SCEPTRE cannot provide (nt_cells at high MOI).
    val control_group

    output:
    // Optional: nothing downstream reads a chunk's MuData (sceptre_chunk_merge takes
    // the tables), and on a screen-scale input each one is gigabytes.
    path "inference_mudata.h5mu", optional: true, emit: inference_mudata
    path "sceptre_per_element_output.tsv.gz", emit: per_element_output
    path "sceptre_per_guide_output.tsv.gz", emit: per_guide_output
    path "sceptre_control_group.json", optional: true, emit: control_group_metadata

    script:
    """
    # Stays at the global cap of 1 BLAS thread: this forks n_processors = task.cpus
    # R workers, and each would otherwise bring its own 64-thread pool.
    cat <<EOF > args.txt
    ${mudata_fp}
    ${params.INFERENCE_SCEPTRE_side}
    ${params.INFERENCE_SCEPTRE_grna_integration_strategy}
    ${params.INFERENCE_SCEPTRE_resampling_approximation}
    ${control_group}
    ${params.INFERENCE_SCEPTRE_resampling_mechanism}
    ${task.cpus}
    EOF

    export SCEPTRE_WRITE_MUDATA=${params.INFERENCE_SCEPTRE_WRITE_MUDATA}
    inference_sceptre.R args.txt
    mv per_element_output.tsv.gz sceptre_per_element_output.tsv.gz
    mv per_guide_output.tsv.gz sceptre_per_guide_output.tsv.gz
    """
}
