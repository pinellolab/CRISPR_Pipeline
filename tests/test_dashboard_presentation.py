import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / 'bin/dashboard_presentation.py'
SPEC = importlib.util.spec_from_file_location('presentation', SCRIPT)
presentation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(presentation)


def report(status):
    return ('<html><head></head><body><div class="live '+status+'">'+status+'</div>'
        '<div class="metric-label">Failed</div><label class="node failed"><span class="node-count">1 failed</span></label>'
        '<section class="failure-evidence"><h2>Failure evidence</h2><details open><summary>mapper</summary>'
        '<pre>original error retained</pre></details></section>'
        '<div class="family-panels"><section id="scientific">Scrublet skipped: threshold unavailable</section>'
        '<figure><img src="data:image/png;base64,AAAA" alt="RNA QC"><figcaption>RNA QC</figcaption></figure></div>'
        '<footer>footer</footer></body></html>')


def test_completed_errors_are_collapsed_separate_history_without_hiding_qc():
    result = presentation.enhance_dashboard_html(report('completed'))
    assert 'id="execution-history"' in result
    assert '<details id="execution-history" class="execution-history">' in result
    assert 'class="failure-evidence"' not in result
    assert result.index('id="execution-history"') > result.index('id="scientific"')
    assert 'original error retained' in result
    assert 'Scrublet skipped: threshold unavailable' in result
    assert 'recovery is unconfirmed' in result
    assert 'Earlier task exceptions' in result
    assert 'class="node review"' in result


def test_active_failure_is_still_prominent():
    result = presentation.enhance_dashboard_html(report('failed'))
    assert 'class="failure-evidence"' in result and '<details open>' in result
    assert 'id="execution-history"' not in result
    assert 'class="node failed"' in result


def test_viewer_is_self_contained_accessible_and_does_not_duplicate_images():
    result = presentation.enhance_dashboard_html(report('completed'))
    assert '<dialog id="plot-viewer"' in result
    assert 'Zoom in' in result and 'Zoom out' in result and 'Actual size' in result
    assert "event.key==='Enter'" in result and 'opener.focus()' in result
    assert result.count('data:image/png;base64,AAAA') == 1
    assert 'https://' not in result
    assert presentation.enhance_dashboard_html(result) == result


def test_unknown_status_is_not_assumed_successful():
    result = presentation.enhance_dashboard_html(report('unknown'))
    assert 'class="failure-evidence"' in result
    assert 'id="execution-history"' not in result


def test_v1_upgrade_keeps_one_viewer_and_is_idempotent():
    old = presentation.enhance_dashboard_html(report('completed')).replace(presentation.MARKER, '<!-- dashboard-presentation-v1 -->')
    updated = presentation.enhance_dashboard_html(old)
    assert updated.count('<dialog id="plot-viewer"') == 1
    assert updated.count('id="execution-history"') == 1
    assert presentation.enhance_dashboard_html(updated) == updated
