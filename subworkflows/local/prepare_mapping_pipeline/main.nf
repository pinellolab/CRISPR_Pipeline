nextflow.enable.dsl=2

include { prepare_covariate } from '../../../modules/local/prepare_covariate'

workflow prepare_mapping_pipeline {
    take:
    ch_samples

    main:
    // Prepare covariate list from the samples channel
    covariate_list = ch_samples
        .map { meta, _reads ->
            [
                modality: meta.modality.toString().toLowerCase(),
                batch: meta.measurement_sets.toString()
            ]
        }
        .unique { row -> "${row.modality}:${row.batch}" }
        // .groupTuple()
        .collect()
        .map { it ->
            def batches_by_modality = [:].withDefault { [] }
            it.each { row ->
                if (!batches_by_modality[row.modality].contains(row.batch)) {
                    batches_by_modality[row.modality] << row.batch
                }
            }

            // Pair modality-specific measurement set IDs by their input order.
            def batch_to_concat_batch = [:]
            batches_by_modality.each { modality, batches ->
                batches.eachWithIndex { batch, idx ->
                    if (!batch_to_concat_batch.containsKey(batch)) {
                        batch_to_concat_batch[batch] = "sample_${idx}"
                    }
                }
            }

            def sorted_batches = batch_to_concat_batch.keySet().toList().sort { a, b ->
                a.toString() <=> b.toString()
            }
            def json = groovy.json.JsonOutput.toJson([
                batch: sorted_batches,
                concat_batch: sorted_batches.collect { batch_to_concat_batch[it] }
            ])
            return json
        }
        .view {"Covariate_list: $it"}

    Prepare_covariate = prepare_covariate(covariate_list)

    emit:
    parsed_covariate_file =  Prepare_covariate.parsed_covariate_file
    covariate_string =  Prepare_covariate.covariate_string
}
