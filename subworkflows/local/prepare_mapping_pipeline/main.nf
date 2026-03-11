nextflow.enable.dsl=2

include { prepare_covariate } from '../../../modules/local/prepare_covariate'

workflow prepare_mapping_pipeline {
    take:
    ch_samples

    main:
    // Prepare covariate list from the samples channel
    covariate_list = ch_samples
        .map { meta, _reads ->
            [meta.measurement_sets, meta.sequencing_run]
        }
        .unique()
        // .groupTuple()
        .map { measurement_sets, sequencing_run ->
            [batch: measurement_sets, cov1: sequencing_run]
        }
        .collect()
        .map { it ->
            def sorted_covariates = it.sort { a, b ->
                a.batch.toString() <=> b.batch.toString() ?: a.cov1.toString() <=> b.cov1.toString()
            }
            def json = groovy.json.JsonOutput.toJson([batch: sorted_covariates.batch.flatten(), cov1: sorted_covariates.cov1.flatten()])
            return json
        }
        .view {"Covariate_list: $it"}

    Prepare_covariate = prepare_covariate(covariate_list)

    emit:
    parsed_covariate_file =  Prepare_covariate.parsed_covariate_file
    covariate_string =  Prepare_covariate.covariate_string
}
