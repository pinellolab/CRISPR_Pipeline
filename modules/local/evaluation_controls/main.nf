
process evaluation_controls {

    cache 'lenient'
    input:

    path mdata


    output:
    path "plots" , emit: evaluation_controls

    script:
            """
            export MPLCONFIGDIR="\$PWD/.matplotlib"
            export XDG_CACHE_HOME="\$PWD/.cache"
            mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME/fontconfig"
            evaluate_controls.py ${mdata}
            """
}
