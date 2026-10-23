# Temporary benchmark dashboards

Five completed report snapshots, published temporarily on **October 23, 2026**
for review. Open a report to browse its QC, guide assignment, inference, and
evaluation sections. These are static snapshots, not live server connections.
They are **not a final manuscript benchmark release**.

<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:16px;margin:24px 0">
<article style="border:1px solid #ccd6e0;border-radius:10px;padding:18px"><h2>Gary Hon</h2><p>Completed calibration dashboard.</p><p><a href="../temporary-dashboards-20261023/GaryHonDataset.html" target="_blank" rel="noopener">Open dashboard ↗</a></p></article>
<article style="border:1px solid #ccd6e0;border-radius:10px;padding:18px"><h2>Engreitz cc-Perturb-seq</h2><p>Completed October 23 dashboard.</p><p><a href="../temporary-dashboards-20261023/Engreitz_ccPerturb.html" target="_blank" rel="noopener">Open dashboard ↗</a></p></article>
<article style="border:1px solid #ccd6e0;border-radius:10px;padding:18px"><h2>HuangFu</h2><p>Completed calibration dashboard, run 04.</p><p><a href="../temporary-dashboards-20261023/HuangFuDataset.html" target="_blank" rel="noopener">Open dashboard ↗</a></p></article>
<article style="border:1px solid #ccd6e0;border-radius:10px;padding:18px"><h2>Gersbach HTV2</h2><p>October 23 rerun with both RNA MAD filters off.</p><p><a href="../temporary-dashboards-20261023/Gersbach_HTV2.html" target="_blank" rel="noopener">Open dashboard ↗</a></p></article>
<article style="border:1px solid #ccd6e0;border-radius:10px;padding:18px"><h2>Gersbach gemX v3</h2><p>October 23 rerun with both RNA MAD filters off.</p><p><a href="../temporary-dashboards-20261023/Gersbach_gemX_v3.html" target="_blank" rel="noopener">Open dashboard ↗</a></p></article>
</div>

Reports are roughly 35–42 MB each because figures are embedded. Open one at a
time; the first load may take a moment. Original server results and Drive copies
are unchanged. These snapshots may predate the newest task-diagnostics feature.

## Interpretation notes

The five runs do not all use identical assay-specific filtering settings.
In particular, the HTV2 and gemX v3 reports are the **both-MAD-filters-off**
reruns, not the earlier MAD-filtered baselines.

Scrublet was enabled in those two reruns but skipped one measurement set in
each because no automatic threshold was resolved. They should not be described
as completely doublet-cleaned datasets. Increased direct-target discrimination
also does not by itself establish better null calibration.

## Temporary hosting and removal

Only the small gallery and deployment instructions are tracked in Git. The
large HTML snapshots come from the temporary GitHub prerelease
`temporary-benchmark-dashboards-20261023`; their compressed archive is checked
against a pinned SHA-256 digest during deployment.

When review is over, remove this page, its navigation/landing-page links, its
entry in `scripts/build_docs_site.py`, the temporary installation step in the
documentation workflow, and `scripts/install_temporary_benchmark_dashboards.py`.
Redeploy Pages, verify the five report URLs are gone, then delete the temporary
prerelease and its assets/tag. The original scientific results must be retained.
Removing public hosting cannot recall already downloaded or cached copies.
