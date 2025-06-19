nextflow.enable.dsl=2

include { prepare_covariate } from '../../../modules/local/prepare_covariate'

workflow prepare_mapping_pipeline {
    take:
    ch_samples

    main:
    // Prepare covariate list from the samples channel
    covariate_list = ch_samples
        .map { meta, _reads ->
            [meta.measurement_sets, meta.lane]
        }
        .unique()
        .groupTuple()
        .map { measurement_sets, lanes ->
            [batch: measurement_sets, cov1: lanes]
        }
        .collect()
        .map { it ->
            def json = groovy.json.JsonOutput.toJson([batch: it.batch.flatten(), cov1: it.cov1.flatten()])
            return json
        }

    Prepare_covariate = prepare_covariate(covariate_list)

    emit:
    parsed_covariate_file =  Prepare_covariate.parsed_covariate_file
    covariate_string =  Prepare_covariate.covariate_string
}
