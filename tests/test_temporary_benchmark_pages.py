import hashlib
import importlib.util
import io
import tarfile
from pathlib import Path

import pytest

PATH = Path(__file__).parents[1] / 'scripts/install_temporary_benchmark_dashboards.py'
SPEC = importlib.util.spec_from_file_location('temporary_pages', PATH)
pages = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pages)


def make_archive(path, names, symlink=False):
    with tarfile.open(path, 'w:gz') as archive:
        for name in names:
            data = b'<!doctype html><html><body>dashboard</body></html>'
            info = tarfile.TarInfo(name)
            if symlink:
                info.type = tarfile.SYMTYPE
                info.linkname = '/etc/passwd'
                archive.addfile(info)
            else:
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))


def test_install_five_reports_and_manifest(tmp_path, monkeypatch):
    archive = tmp_path / 'reports.tar.gz'
    make_archive(archive, sorted(pages.FILES))
    monkeypatch.setattr(pages, 'DIGEST', hashlib.sha256(archive.read_bytes()).hexdigest())
    pages.install(archive, tmp_path / 'site')
    target = tmp_path / 'site/temporary-dashboards-20261023'
    assert {p.name for p in target.iterdir()} == pages.FILES | {'manifest.json'}


@pytest.mark.parametrize('symlink', [False, True])
def test_reject_paths_and_symlinks(tmp_path, monkeypatch, symlink):
    archive = tmp_path / 'reports.tar.gz'
    names = sorted(pages.FILES) if symlink else [f'../{n}' for n in sorted(pages.FILES)]
    make_archive(archive, names, symlink=symlink)
    monkeypatch.setattr(pages, 'DIGEST', hashlib.sha256(archive.read_bytes()).hexdigest())
    with pytest.raises(ValueError):
        pages.install(archive, tmp_path / 'site')
    assert not (tmp_path / 'site').exists()


def test_reject_checksum_mismatch(tmp_path):
    archive = tmp_path / 'reports.tar.gz'
    make_archive(archive, sorted(pages.FILES))
    with pytest.raises(ValueError, match='checksum'):
        pages.install(archive, tmp_path / 'site')
