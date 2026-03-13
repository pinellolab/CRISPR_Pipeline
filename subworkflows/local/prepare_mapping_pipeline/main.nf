nextflow.enable.dsl=2

include { prepare_covariate } from '../../../modules/local/prepare_covariate'

workflow prepare_mapping_pipeline {
    take:
    ch_samples

    main:
    // Prepare covariate list from the samples channel
    covariate_list = ch_samples
        .map { meta, _reads ->
            [meta.measurement_sets]
        }
        .unique()
        // .groupTuple()
        .map { measurement_sets ->
            [batch: measurement_sets]
        }
        .collect()
        .map { it ->
            def sorted_covariates = it.sort { a, b ->
                a.batch.toString() <=> b.batch.toString()
            }
            def json = groovy.json.JsonOutput.toJson([batch: sorted_covariates.batch.flatten()])
            return json
        }
        .view {"Covariate_list: $it"}

    Prepare_covariate = prepare_covariate(covariate_list)

    emit:
    parsed_covariate_file =  Prepare_covariate.parsed_covariate_file
    covariate_string =  Prepare_covariate.covariate_string
}
