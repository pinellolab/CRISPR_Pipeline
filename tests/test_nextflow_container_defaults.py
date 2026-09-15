import re
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
    assert "crispr_pipeline/perturbo:v2-cis-mask" not in config


def test_perturbo_is_pinned_by_digest():
    """:v2-dev is republished on every push to v2-port, so a tag pin can change
    under a running test. Assert the invariant, not a literal that goes stale."""
    config = (REPO_ROOT / "nextflow.config").read_text()
    pin = re.search(r"perturbo\s*=\s*'([^']+)'", config).group(1)

    assert pin.startswith("ghcr.io/pinellolab/perturbo@sha256:"), pin
    assert len(pin.split("sha256:")[1]) == 64, pin


def test_benchmark_assets_are_resolved_from_the_pipeline_checkout():
    config = (REPO_ROOT / "nextflow.config").read_text()

    assert 'ENCODE_BED_DIR = "${projectDir}/encode_bed_files"' in config


def test_library_thread_pools_are_capped_inside_tasks():
    """Every task must cap OpenBLAS/OpenMP/Arrow threads, or pools size to the host.

    The `env` scope is the vector that reaches inside the container on the local
    and SLURM executors; `beforeScript` runs on the host and is wiped before the
    container starts. A process that does real parallel work opts back in with
    task.cpus in its own script block.
    """
    # nextflow_cc.config stands alone -- it never includes nextflow.config -- so the
    # CC-Perturb-seq path only has the caps if this file carries its own copy.
    for name in ("nextflow.config", "nextflow_cc.config"):
        config = (REPO_ROOT / name).read_text()
        env_block = re.search(r"^env \{(.*?)^\}", config, re.S | re.M)
        assert env_block, f"{name} has no top-level env {{}} scope"
        body = env_block.group(1)
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_MAX_THREADS", "POLARS_MAX_THREADS"):
            assert re.search(rf"{var}\s*=\s*'1'", body), f"{var} not capped to 1 in {name}"

    assignment = (REPO_ROOT / "modules/local/guide_assignment_sceptre/main.nf").read_text()
    assert "export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus}" in assignment
