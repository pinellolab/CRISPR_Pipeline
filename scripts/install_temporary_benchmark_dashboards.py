#!/usr/bin/env python3
"""Install only the five approved temporary Pages reports, never arbitrary tar paths."""

import argparse
import hashlib
import json
import tarfile
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from dashboard_presentation import enhance_dashboard_html

DIGEST = "cecd79bc7a072cfe1eb97c934b3391aef375156184cd0bd50005f5b0b22a3475"
FILES = {
    "Engreitz_ccPerturb.html", "GaryHonDataset.html", "HuangFuDataset.html",
    "Gersbach_HTV2.html", "Gersbach_gemX_v3.html",
}


def install(archive_path, site_dir):
    digest = hashlib.sha256()
    with archive_path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest() != DIGEST:
        raise ValueError('Temporary dashboard archive checksum mismatch')
    target = site_dir / 'temporary-dashboards-20261023'
    with tarfile.open(archive_path, 'r:gz') as archive:
        members = archive.getmembers()
        if len(members) != 5 or {m.name for m in members} != FILES:
            raise ValueError('Archive must contain exactly the five approved reports')
        if any(not m.isfile() or m.size > 100_000_000 for m in members):
            raise ValueError('Only regular HTML files under 100 MB are allowed')
        target.mkdir(parents=True, exist_ok=True)
        manifest = []
        for member in members:
            payload = archive.extractfile(member).read()
            if b'<html' not in payload[:1024].lower():
                raise ValueError(f'Not an HTML dashboard: {member.name}')
            payload = enhance_dashboard_html(payload.decode('utf-8')).encode('utf-8')
            (target / member.name).write_bytes(payload)
            manifest.append(dict(file=member.name, bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest()))
        (target / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Installed {len(manifest)} temporary dashboards in {target}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--site-dir', type=Path, default=Path('build/site'))
    args = parser.parse_args()
    install(args.archive, args.site_dir)
