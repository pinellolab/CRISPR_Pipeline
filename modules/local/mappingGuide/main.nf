
process mappingGuide {
    cache 'lenient'
    debug true

    input:
    tuple val(meta), path(reads)
    path guide_index
    path t2g_guide
    path parsed_seqSpec_file
    path barcode_file

    output:
    path "*_ks_guide_out", emit: ks_guide_out_dir
    path "*_ks_guide_out/counts_unfiltered/adata.h5ad", emit: ks_guide_out_adata

    script:
        def batch = meta.measurement_sets
        def fastq_files = reads.join(' ')
        """
        echo "Processing $batch with $fastq_files"

        k_bin=\$(type -p kallisto)
        bustools_bin=\$(type -p bustools)
        chemistry=\$(extract_parsed_seqspec.py --file ${parsed_seqSpec_file})

        kb count -i ${guide_index} -g ${t2g_guide} --verbose -w ${barcode_file} \\
                --h5ad --kallisto \$k_bin --bustools \$bustools_bin -x \$chemistry -o ${batch}_ks_guide_out -t ${task.cpus} \\
                ${fastq_files} --overwrite

        echo "gRNA KB mapping Complete"
        """
}
