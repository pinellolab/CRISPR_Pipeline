from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_every_process_has_a_base_container_fallback():
    config = (REPO_ROOT / "nextflow.config").read_text()
    process_start = config.index("// Global process configuration")
    first_override = config.index("withName:", process_start)
    global_process_config = config[process_start:first_override]

    assert "container = { params.containers.base }" in global_process_config


def test_specialized_processes_retain_their_container_overrides():
    config = (REPO_ROOT / "nextflow.config").read_text()

    assert "container = { params.containers.bediting }" in config
    assert "container = { params.containers.cleanser }" in config
    assert "container = { params.containers.sceptre }" in config
    assert "container = { params.containers.perturbo }" in config


def test_benchmark_assets_are_resolved_from_the_pipeline_checkout():
    config = (REPO_ROOT / "nextflow.config").read_text()

    assert 'ENCODE_BED_DIR = "${projectDir}/encode_bed_files"' in config
