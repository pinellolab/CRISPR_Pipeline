process mappingGuideBaseEditingFlash {
    cache 'lenient'
    debug true
    stageOutMode 'copy'

    input:
    tuple val(meta), path(reads)
    path parsed_seqSpec_file
    path barcode_file
    val  is_10xv3v
    path metadata_file
    val reverse_complement_guides
    val spacer_tag

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
        mkdir -p ${batch}_ks_guide_out

        flash_base_editing_mapping.py \\
        --barcode_inclusion_list_fn ${barcode_file} \\
        --fastq "${fastq_files}" \\
        --guide_set_fn ${metadata_file} \\
        --chemistry "\$chemistry" \\
        --cores ${task.cpus} \\
        --downsample_reads 0 \\
        --output_prefix ${batch}_ks_guide_out \\
        --reverse_complement_guides ${reverse_complement_guides} \\
        --spacer_tag "${spacer_tag}" \\
        --tolerance ${params.BASEEDITING_FLASH_tolerance} \\
        --guide_len ${params.BASEEDITING_FLASH_guide_len} \\
        --device ${params.BASEEDITING_FLASH_device} \\
        --fastq_chunk_size ${params.BASEEDITING_FLASH_fastq_chunk_size} \\
        --gpu_read_chunk ${params.BASEEDITING_FLASH_gpu_read_chunk}
        """
}
