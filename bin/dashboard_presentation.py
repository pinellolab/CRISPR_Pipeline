#!/usr/bin/env python3
"""Presentation-only upgrades for new and already-completed HTML dashboards."""

import argparse
import re
from pathlib import Path

MARKER = '<!-- dashboard-presentation-v1 -->'

VIEWER = r'''
<style id="dashboard-presentation-style">
figure img { cursor: zoom-in; }
.execution-history { margin: 18px 0; padding: 16px; border: 1px solid #cbd5e1; border-radius: 12px; background: #f8fafc; }
.execution-history>summary { color:#334155; font-size:17px; }
.execution-history details { border-color:#cbd5e1; }
.execution-history summary { color:#334155; }
.execution-history p { color:#475569; }
.run-completed .node.review .node-status,.run-completed .review.status-badge { background:#64748b; }
.run-completed .pill.failed,.run-completed .pill.aborted,.run-completed .pill.interrupted { background:#f1f5f9; color:#475569; }
#plot-viewer { border:0; border-radius:12px; padding:12px; width:96vw; max-width:96vw; height:94vh; max-height:94vh; background:#fff; color:#0f172a; }
#plot-viewer::backdrop { background:rgba(15,23,42,.8); }
.plot-toolbar { display:flex; align-items:center; flex-wrap:wrap; gap:8px; padding-bottom:10px; }
.plot-toolbar button,.plot-toolbar a { padding:8px 12px; border:1px solid #94a3b8; border-radius:6px; color:#0f172a; background:#fff; cursor:pointer; font:14px system-ui; }
#plot-viewer-title { flex:1; min-width:140px; overflow-wrap:anywhere; }
#plot-viewport { width:100%; height:calc(100% - 60px); overflow:auto; background:#eef2f6; }
#plot-canvas { display:flex; align-items:flex-start; justify-content:center; min-width:100%; min-height:100%; }
#plot-viewer-image { max-width:none; flex:none; display:block; margin:0; cursor:default; }
.plot-help { font-size:12px; color:#475569; }
</style>
<dialog id="plot-viewer" aria-labelledby="plot-viewer-title">
 <div class="plot-toolbar">
  <strong id="plot-viewer-title">Plot</strong>
  <button type="button" data-plot-action="out" aria-label="Zoom out">−</button>
  <output id="plot-zoom" aria-live="polite">100%</output>
  <button type="button" data-plot-action="in" aria-label="Zoom in">+</button>
  <button type="button" data-plot-action="fit">Fit</button>
  <button type="button" data-plot-action="actual">Actual size</button>
  <a id="plot-download" download="plot.png">Download</a>
  <button type="button" data-plot-action="close">Close</button>
 </div>
 <div id="plot-viewport" tabindex="0" aria-label="Zoomed plot; scroll to inspect">
  <div id="plot-canvas"><img id="plot-viewer-image" alt=""></div>
 </div>
</dialog>
<script id="dashboard-plot-viewer-script">
(() => {
 const dialog=document.getElementById('plot-viewer');
 const image=document.getElementById('plot-viewer-image');
 const viewport=document.getElementById('plot-viewport');
 const canvas=document.getElementById('plot-canvas');
 const label=document.getElementById('plot-zoom');
 let scale=1, opener=null;
 function zoom(value) {
   scale=Math.max(.05,Math.min(16,value));
   const width=Math.max(1,image.naturalWidth*scale);
   const height=Math.max(1,image.naturalHeight*scale);
   image.style.width=width+'px'; image.style.height=height+'px';
   canvas.style.width=Math.max(viewport.clientWidth,width)+'px';
   canvas.style.height=Math.max(viewport.clientHeight,height)+'px';
   label.value=Math.round(scale*100)+'%';
 }
 function fit() { zoom(Math.min(1,(viewport.clientWidth-16)/image.naturalWidth,(viewport.clientHeight-16)/image.naturalHeight)); viewport.scrollTo(0,0); }
 function openPlot(plot) {
   opener=plot;
   document.getElementById('plot-viewer-title').textContent=plot.alt || 'Plot';
   image.alt=plot.alt || 'Expanded plot';
   image.onload=fit;
   image.src=plot.currentSrc || plot.src;
   document.getElementById('plot-download').href=image.src;
   if (!dialog.open) dialog.showModal();
   if (image.complete && image.naturalWidth) fit();
 }
 document.querySelectorAll('figure img').forEach(plot => {
   plot.tabIndex=0; plot.setAttribute('role','button');
   plot.setAttribute('aria-label','Expand plot: '+(plot.alt || 'plot'));
   plot.title='Click to enlarge and zoom';
   plot.addEventListener('click',() => openPlot(plot));
   plot.addEventListener('keydown',event => { if (event.key==='Enter' || event.key===' ') { event.preventDefault(); openPlot(plot); } });
 });
 dialog.addEventListener('click',event => {
   const action=event.target.dataset.plotAction;
   if(action==='in') zoom(scale*1.25);
   if(action==='out') zoom(scale/1.25);
   if(action==='fit') fit();
   if(action==='actual') zoom(1);
   if(action==='close') dialog.close();
 });
 dialog.addEventListener('close',() => { if(opener) opener.focus(); });
 viewport.addEventListener('wheel',event => { if(event.ctrlKey || event.metaKey) { event.preventDefault(); zoom(scale*(event.deltaY<0?1.15:1/1.15)); } },{passive:false});
 dialog.addEventListener('keydown',event => { if(event.key==='+' || event.key==='=') {event.preventDefault();zoom(scale*1.25);} if(event.key==='-') {event.preventDefault();zoom(scale/1.25);} });
 window.addEventListener('resize',() => { if(dialog.open) fit(); });
})();
</script>
'''


def enhance_dashboard_html(document, run_status=None):
    """Retain all evidence and scientific warnings; change only its presentation."""
    if MARKER in document:
        return document
    if run_status is None:
        match = re.search(r'class="live\s+([^" ]+)', document)
        run_status = match.group(1) if match else 'unknown'
    completed = run_status.lower() == 'completed'
    if completed:
        evidence = []
        def move_evidence(match):
            block = match.group(0)
            block = re.sub(r'<details\s+open>', '<details>', block)
            block = block.replace('class="failure-evidence"', 'class="historical-evidence"')
            block = block.replace('<h2>Failure evidence</h2>', '<h3>Earlier execution diagnostics</h3>')
            evidence.append(block)
            return ''
        document = re.sub(r'<section class="failure-evidence">.*?</section>', move_evidence, document, flags=re.S)
        # A successful workflow may have ignored failures, so never claim that
        # every historical exception recovered without task-level evidence.
        if evidence:
            history = ('<details id="execution-history" class="execution-history">'
                '<summary>Execution history &amp; diagnostics — completed workflow</summary>'
                '<p>This workflow completed successfully. Earlier errors and interruptions are retained here, '
                'not presented as an active pipeline failure. A task is recovered only where a later '
                'successful attempt is documented; otherwise its recovery is unconfirmed. '
                'Result-quality and skipped-analysis warnings remain in their scientific sections.</p>'
                + ''.join(evidence) + '</details>')
            document = document.replace('<footer>', history+'<footer>', 1)
        document = document.replace('<body>', '<body class="run-completed">', 1)
        document = document.replace('class="node failed"', 'class="node review"')
        document = document.replace('class="status-badge failed">failed', 'class="status-badge review">task exceptions')
        document = document.replace('<div class="metric-label">Failed</div>', '<div class="metric-label">Earlier task exceptions</div>')
        document = re.sub(r'(<span class="node-count">[^<]*) failed', r'\1 earlier exceptions', document)
    return document.replace('</body>', MARKER+VIEWER+'</body>', 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error('Use a separate output file to preserve the original report')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(enhance_dashboard_html(args.input.read_text(encoding='utf-8')), encoding='utf-8')
