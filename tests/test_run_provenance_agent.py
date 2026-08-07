import json
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "run_provenance_agent.py"


def run(*args, check=True):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        text=True,
        capture_output=True,
        check=check,
    )


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def make_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    git(repo, "init")
    git(repo, "config", "user.email", "test@example.org")
    git(repo, "config", "user.name", "Test User")
    git(repo, "remote", "add", "origin", str(repo))
    (repo / "main.nf").write_text("workflow { log.info 'test' }\n")
    (repo / "nextflow.config").write_text("params.outdir = 'results'\n")
    git(repo, "add", "main.nf", "nextflow.config")
    git(repo, "commit", "-m", "fixture")
    source = tmp_path / "source"
    shutil.copytree(repo, source, ignore=shutil.ignore_patterns(".git"))
    return repo, source


def test_prepare_and_check_clean_source(tmp_path):
    repo, source = make_repo(tmp_path)
    artifact = tmp_path / "params.json"
    artifact.write_text('{"input":"sample.csv"}\n')
    output = tmp_path / "results" / "pipeline_info"

    run(
        "prepare",
        "--pipeline-repo",
        repo,
        "--source-dir",
        source,
        "--output-dir",
        output,
        "--run-name",
        "fixture",
        "--artifact",
        f"params={artifact}",
    )
    data = json.loads((output / "repository_provenance.json").read_text())

    assert data["pipeline_source"]["source_matches_expected_commit"] is True
    assert data["pipeline_source"]["expected_commit"] == git(repo, "rev-parse", "HEAD")
    assert data["launch_artifacts"]["params"]["sha256"]
    assert "manifest {" in (output / "provenance.generated.config").read_text()
    run(
        "check",
        "--pipeline-repo",
        repo,
        "--metadata",
        output / "repository_provenance.json",
    )


def test_check_fails_when_recorded_config_changes(tmp_path):
    repo, source = make_repo(tmp_path)
    config = tmp_path / "run.config"
    config.write_text("process.maxForks = 4\n")
    output = tmp_path / "pipeline_info"
    run(
        "prepare",
        "--pipeline-repo",
        repo,
        "--source-dir",
        source,
        "--output-dir",
        output,
        "--run-name",
        "fixture",
        "--artifact",
        f"config={config}",
    )
    config.write_text("process.maxForks = 8\n")

    result = run(
        "check",
        "--pipeline-repo",
        repo,
        "--metadata",
        output / "repository_provenance.json",
        check=False,
    )

    assert result.returncode == 2
    assert "artifact config changed" in result.stderr


def test_prepare_rejects_source_that_does_not_match_commit(tmp_path):
    repo, source = make_repo(tmp_path)
    (source / "main.nf").write_text("workflow { error 'changed' }\n")

    result = run(
        "prepare",
        "--pipeline-repo",
        repo,
        "--source-dir",
        source,
        "--output-dir",
        tmp_path / "pipeline_info",
        "--run-name",
        "fixture",
        check=False,
    )

    assert result.returncode == 2
    assert "executing source does not match commit" in result.stderr
