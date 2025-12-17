process mappingGuideBaseEditing {
    cache 'lenient'
    debug true
    stageOutMode 'copy'

    input:
    tuple val(meta), path(reads)
    path parsed_seqSpec_file
    path barcode_file
    val  is_10xv3v   // should be "true" or "false"
    path metadata_file
    

    output:
    path "*_ks_guide_out", emit: ks_guide_out_dir
    path "*_ks_guide_out/counts_unfiltered/adata.h5ad", emit: ks_guide_out_adata

    script:
        def batch       = meta.measurement_sets
        def fastq_files = reads.join(' ')
        """
        set -euo pipefail
        echo "Processing $batch with $fastq_files"

        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})
        echo "chemistry"

        mkdir -p ${batch}_ks_guide_out
        #cd ${batch}_ks_guide_out
        echo "base editing script will start"
        pip list

        base_editing_mapping.py \\
        --barcode_inclusion_list_fn ${barcode_file} \\
        --fastq "${fastq_files}" \\
        --guide_set_fn ${metadata_file} \\
        --chemistry "\$chemistry" \\
        --cores ${task.cpus} \\
        --downsample_reads 0 \\
        --output_prefix ${batch}_ks_guide_out
        echo "gRNA KB mapping Complete"
        """
}
