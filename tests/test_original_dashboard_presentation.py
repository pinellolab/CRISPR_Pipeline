"""Exercise the actual original generators, not just the retrofit helper."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'bin'))


@pytest.mark.parametrize('script', ['create_dashboard.py', 'create_dashboard_HASHING.py'])
def test_original_generators_embed_full_plot_viewer(script, monkeypatch):
    monkeypatch.chdir(ROOT)
    spec = importlib.util.spec_from_file_location(script[:-3], ROOT / 'bin' / script)
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    frame = pd.DataFrame([dict(modality='scRNA', subject='RNA QC', description='QC figure',
        value_display='Scientific metric: 42', highlighted=False, table=pd.DataFrame(),
        table_description='', image=['figures/example.png'], image_description=['Original plot'])])
    result = generator.df_to_html_side_tabs(frame, {})
    assert 'Error generating dashboard' not in result
    assert 'Scientific metric: 42' in result
    assert 'data-imgsrc="figures/example.png"' in result
    assert '<dialog id="plot-viewer"' in result
    assert "closest('button[data-imgsrc]')" in result
    assert 'event.stopImmediatePropagation()' in result
    assert 'id="execution-history"' not in result  # no invented workflow status
