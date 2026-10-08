import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('qc_dashboard', ROOT/'bin/render_wandb_pipeline_dashboard.py')
dashboard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dashboard)


def test_postconcat_node_order_and_process_routing():
    families = [f[0] for f in dashboard.FAMILIES]
    assert families.index('guide_assignment') + 1 == families.index('postconcat_qc')
    assert families.index('postconcat_qc') + 1 == families.index('inference')
    for name in ('embedding_before_clone', 'embedding_after_clone', 'remove_clonal_cells', 'filter_hto_post_clone'):
        assert dashboard.family_for('pipeline:'+name) == 'postconcat_qc'
    assert dashboard.image_family(Path('pipeline_dashboard/additional_qc/embeddings/part01/embedding_qc/umap_qc_panel.png')) == 'postconcat_qc'


def test_same_plot_kept_once_per_stage_and_omissions_reported(tmp_path):
    for stage in ('before_clone', 'after_clone'):
        path = tmp_path/stage
        path.mkdir()
        (path/'umap_qc_panel.png').write_bytes(b'identical-png')
        (path/'embedding_qc_metrics.json').write_text(json.dumps({'stage': stage}))
    selected = dashboard.collect_images(tmp_path)
    assert len(selected['postconcat_qc']) == 2
    omitted = []
    dashboard.collect_images(tmp_path, 1, omitted)
    assert len(omitted) == 2
    gallery = dashboard.image_gallery(selected['postconcat_qc'])
    assert 'before_clone umap_qc_panel' in gallery
    assert 'after_clone umap_qc_panel' in gallery


def test_stage_metrics_are_read_without_final_dashboard(tmp_path):
    path = tmp_path/'postconcat_embedding_qc'/'before_clone'/'embedding_qc'
    path.mkdir(parents=True)
    (path/'embedding_qc_metrics.json').write_text(json.dumps({'stage':'before_clone', 'cells':43, 'effective_pcs':12}))
    content = dashboard.postconcat_qc_content(tmp_path)
    assert 'before_clone' in content and '43' in content and '12' in content


def test_live_publication_matches_emitted_directories():
    for module, folder in (('postconcat_embedding_qc','embedding_qc'), ('remove_clonal_cells','clone_qc'),
                           ('filter_guide_assignment_qc','guide_assignment_qc'), ('filter_hto_post_clone','hto_qc')):
        text = (ROOT/'modules/local'/module/'main.nf').read_text()
        assert f"pattern: '{folder}'" in text or f'pattern: "{folder}"' in text
        assert f'{folder}/*' not in text
