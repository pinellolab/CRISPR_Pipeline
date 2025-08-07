process downloadReference {
    cache 'lenient'

    input:
    val ref_name
    val use_igvf_reference

    output:
    path "transcriptome_index.idx", emit: transcriptome_idx
    path "transcriptome_t2g.txt", emit: t2g_transcriptome_index

    script:
    def igvf_url = "https://api.data.igvf.org/reference-files/IGVFFI9561BASO/@@download/IGVFFI9561BASO.tar.gz"

    if (use_igvf_reference) {
        """
        wget -qO igvf_reference.tar.gz ${igvf_url}
        mkdir -p igvf_extracted
        tar -xzf igvf_reference.tar.gz -C igvf_extracted
        mv igvf_extracted/*.idx transcriptome_index.idx
        mv igvf_extracted/t2g.txt transcriptome_t2g.txt
        """
    } else {
        """
        k_bin=\$(type -p kallisto)
        kb ref -d ${ref_name} -i transcriptome_index.idx -g transcriptome_t2g.txt --kallisto \$k_bin
        """
    }
}
