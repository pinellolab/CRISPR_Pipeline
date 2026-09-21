// Tiny containerized publication regression test; no sequencing data.
nextflow.enable.dsl=2
params.outdir = 'publication_smoke'

process emit_qc {
    publishDir "${params.outdir}/live", mode: 'copy', pattern: 'embedding_qc'
    output:
    path 'embedding_qc', emit: qc
    script:
    """
    mkdir -p embedding_qc
    echo QC_PUBLICATION_SENTINEL > embedding_qc/test_plot.txt
    """
}

process collect_qc {
    publishDir "${params.outdir}/collected", mode: 'copy'
    input:
    path qc
    output:
    path 'collected_qc'
    script:
    """
    test -L ${qc}
    find -L ${qc} -mindepth 1 -type f -print -quit | grep -q .
    mkdir collected_qc
    cp -R ${qc}/. collected_qc/
    """
}

workflow {
    emit_qc()
    collect_qc(emit_qc.out.qc)
}
