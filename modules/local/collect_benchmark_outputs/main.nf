// Assemble the benchmark_output directory the dashboard copies: the TF
// benchmark's files when the benchmark ran (the placeholder holds only a
// .gitkeep) and the non-targeting calibration check under
// benchmark_output/calibration/. The dashboard ignores files it does not know,
// so nothing there changes.
process collect_benchmark_outputs {
    cache 'lenient'

    input:
        path benchmark_dir, stageAs: 'benchmark_in'
        path calibration_dir, stageAs: 'calibration_in'

    output:
        path 'benchmark_output', emit: collected

    script:
        """
        mkdir -p benchmark_output
        if find -L benchmark_in -mindepth 1 -type f ! -name '.gitkeep' -print -quit | grep -q .; then
            cp -RL benchmark_in/. benchmark_output/
        fi
        if find -L calibration_in -mindepth 1 -type f ! -name '.gitkeep' -print -quit | grep -q .; then
            mkdir -p benchmark_output/calibration
            cp -RL calibration_in/. benchmark_output/calibration/
        fi
        """
}
