"""Keep documentation/source compatibility checks in the normal test suite."""

import subprocess
import sys
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_documentation_and_plugin_contract():
    subprocess.run([sys.executable, str(ROOT / "scripts/check_docs_sync.py")], cwd=ROOT, check=True)


def test_bundled_plugin_installs_idempotently(tmp_path):
    env = dict(os.environ, NXF_PLUGINS_DIR=str(tmp_path / "nextflow-plugins"))
    installer = ROOT / "run_support/install_dashboard_plugin.sh"
    first = subprocess.run(["bash", str(installer)], cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    second = subprocess.run(["bash", str(installer)], cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    assert "Installed dashboard plugin" in first.stdout
    assert "already installed" in second.stdout
    assert (tmp_path / "nextflow-plugins/nf-crispr-dashboard-0.1.0/.archive.sha256").is_file()
