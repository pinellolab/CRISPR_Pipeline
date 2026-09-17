include { prepare_assignment } from '../../../modules/local/prepare_assignment'
include { guide_assignment_cleanser } from '../../../modules/local/guide_assignment_cleanser'
include { guide_assignment_sceptre } from '../../../modules/local/guide_assignment_sceptre'
include { mudata_concat } from '../../../modules/local/mudata_concat'
include { filter_guide_assignment_qc } from '../../../modules/local/filter_guide_assignment_qc'

workflow guide_assignment_pipeline {

    take:
    mudata_input

    main:
    Prepare_assignment = prepare_assignment(mudata_input)

    if (params.GUIDE_ASSIGNMENT_method == "cleanser") {
        Guide_Assignment = guide_assignment_cleanser(
            Prepare_assignment.prepare_assignment_mudata.flatten(), 
            params.GUIDE_ASSIGNMENT_cleanser_probability_threshold, 
            params.GUIDE_ASSIGNMENT_capture_method
        )
        guide_assignment_collected = Guide_Assignment.guide_assignment_mudata_output
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }
        Mudata_concat = mudata_concat(guide_assignment_collected, params.DUAL_GUIDE)
    }
    else if (params.GUIDE_ASSIGNMENT_method == "sceptre") {
        Guide_Assignment = guide_assignment_sceptre(
            Prepare_assignment.prepare_assignment_mudata.flatten(), 
            params.GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold, 
            params.GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep
        )
        guide_assignment_collected = Guide_Assignment.guide_assignment_mudata_output
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }
        Mudata_concat = mudata_concat(guide_assignment_collected, params.DUAL_GUIDE)
    } else {
        error("Invalid GUIDE_ASSIGNMENT_method: ${params.GUIDE_ASSIGNMENT_method}")
    }

    GuideAssignmentQC = filter_guide_assignment_qc(Mudata_concat.concat_mudata)

    emit:
    concat_mudata = GuideAssignmentQC.filtered_mudata
    guide_assignment_qc = GuideAssignmentQC.guide_assignment_qc
}
