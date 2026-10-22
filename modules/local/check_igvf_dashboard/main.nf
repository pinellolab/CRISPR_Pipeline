process check_igvf_dashboard {
    tag 'plugin dashboard required for IGVF package'

    input:
    path legacy_archive

    output:
    path 'igvf_dashboard_present.txt', emit: marker

    script:
    """
    test -s '${params.outdir}/pipeline_dashboard/execution_dashboard.html' || { echo 'IGVF package requires the live plugin execution_dashboard.html' >&2; exit 1; }
    echo 'Plugin execution dashboard present' > igvf_dashboard_present.txt
    """
}
