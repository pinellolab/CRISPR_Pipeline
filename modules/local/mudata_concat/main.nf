// The pipeline's single cell filter lives here, in mudata_concat.py: this step runs
// right after guide assignment, so the cell population it decides is the one both
// inference methods, every derived covariate and both CRT pools inherit.
// require_assigned_guide drops cells with no assigned guide; the pre-filter counts
// are recorded in .uns so the guide-assignment rate stays reportable. Under
// DUAL_GUIDE, collapse_guides then applies its own stricter rule (exactly two
// guides on one element), which subsumes this one rather than repeating it.
process mudata_concat {
    cache 'lenient'

    input:
        path (mudata_input, name: "?/*")
        val fraction_cells
        val require_assigned_guide
        val dual_guide
    output:
        path "concat_mudata.h5mu", emit: concat_mudata

    script:
    """
        python ${projectDir}/bin/mudata_concat.py -i ${mudata_input} -o concat_mudata.h5mu -g ${fraction_cells} --require-assigned-guide ${require_assigned_guide}

        if [ "${dual_guide}" = "true" ]; then
            echo "Dual guide mode enabled, processing accordingly."
            python ${projectDir}/bin/collapse_guides.py concat_mudata.h5mu concat_mudata_collapsed.h5mu
            mv concat_mudata_collapsed.h5mu concat_mudata.h5mu

        else
            echo "Single guide mode, proceeding with standard processing."
        fi
    """
}










