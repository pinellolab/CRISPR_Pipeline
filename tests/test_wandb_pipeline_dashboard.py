import argparse
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "render_wandb_pipeline_dashboard.py"
SPEC = importlib.util.spec_from_file_location("wandb_pipeline_dashboard", SCRIPT)
dashboard = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(dashboard)


def test_render_builds_clickable_family_dashboard(tmp_path):
    trace = tmp_path / "trace.tsv"
    trace.write_text(
        "process\tname\tstatus\trealtime\tpeak_rss\n"
        "NFCORE_CRISPR:CRISPR_PIPELINE:seqspeccheck\tseqspeccheck (sample)\tCOMPLETED\t4s\t20 MB\n"
        "NFCORE_CRISPR:CRISPR_PIPELINE:inference_pipeline:inference_sceptre\tsceptre (1)\tFAILED\t2m\t1 GB\n",
        encoding="utf-8",
    )
    guide_report = tmp_path / "guides.json"
    guide_report.write_text(
        '{"row_count": 42, "targeting_rows": 40, "control_rows": 2}',
        encoding="utf-8",
    )
    seqspec_table = tmp_path / "seqspec.csv"
    seqspec_table.write_text(
        "Sample,Config,TotalHits,HitRatio,PosPurity,FlankPurity,FinalScore,IsWinner\n"
        "sample,config_a,100,0.9,0.8,0.7,12.5,true\n",
        encoding="utf-8",
    )
    seqspec_image = tmp_path / "seqspec.png"
    seqspec_image.write_bytes(b"\x89PNG\r\n\x1a\n")

    result = dashboard.render(
        argparse.Namespace(
            trace=trace,
            run_id="run-123",
            run_name="TAP-seq chr8",
            status="interrupted",
            guide_report=guide_report,
            seqspec_table=seqspec_table,
            seqspec_image=seqspec_image,
        )
    )

    assert "Live dependency view" in result
    assert "selectFamily('seqspec')" in result
    assert 'id="family-inference"' in result
    assert "FAILED" in result
    assert "42" in result
    assert "data:image/png;base64," in result
    assert "No credentials, FASTQs or unbounded task logs embedded" in result
    assert "location.hash.slice(1)" in result


def test_processes_are_mapped_to_expected_families():
    assert dashboard.family_for("pipeline:seqspeccheck") == "seqspec"
    assert dashboard.family_for("pipeline:mapping_rna_pipeline:mappingrna") == "mapping"
    assert dashboard.family_for("pipeline:guide_assignment_sceptre") == "guide_assignment"
    assert dashboard.family_for("pipeline:inference_pipeline:inference_perturbo") == "inference"
    assert dashboard.family_for("pipeline:additional_qc") == "evaluation"
    assert dashboard.family_for("pipeline:dashboard") == "final"


def test_qc_catalog_and_sanitized_failure_evidence_are_embedded(tmp_path):
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / ".exitcode").write_text("1\n", encoding="utf-8")
    (workdir / ".command.err").write_text(
        "starting tool\nAPI_TOKEN=do-not-publish-this-value\nactual failure\n",
        encoding="utf-8",
    )
    trace = tmp_path / "trace.tsv"
    trace.write_text(
        "process\tname\tstatus\texit\tworkdir\n"
        f"pipeline:inference_sceptre\tsceptre (1)\tFAILED\t1\t{workdir}\n",
        encoding="utf-8",
    )
    nextflow_log = tmp_path / "nextflow.log"
    nextflow_log.write_text("ERROR ~ Error executing process > sceptre\n", encoding="utf-8")
    guide_report = tmp_path / "guides.json"
    guide_report.write_text("{}", encoding="utf-8")
    seqspec_table = tmp_path / "seqspec.csv"
    seqspec_table.write_text("Sample,IsWinner\n", encoding="utf-8")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        '{"dashboard_summary":{"filtering_summary":{"final_mudata":{"cells":123,"gene_features":45}}},'
        '"observed_metrics":{"additional_qc":{"guide":{"rows":[{"batch":"all",'
        '"n_cells_with_guide":100,"frac_cells_with_guide":0.8}]}}}}',
        encoding="utf-8",
    )

    result = dashboard.render(
        argparse.Namespace(
            trace=trace, run_id="failed-run", run_name="failed run", status="failed",
            guide_report=guide_report, seqspec_table=seqspec_table, seqspec_image=None,
            qc_metrics_json=metrics, artifact_dir=None, max_image_bytes=1_000_000,
            nextflow_log=nextflow_log, tail_lines=10,
        )
    )

    assert "Failure evidence" in result
    assert "Error executing process" in result
    assert "actual failure" in result
    assert "do-not-publish-this-value" not in result
    assert "API_TOKEN=&lt;redacted&gt;" in result
    assert "Final cells" in result and "123" in result
    assert "Assignment rate" in result and "0.800" in result


def test_final_dashboard_inference_tables_are_bounded_and_searchable(tmp_path):
    final_dashboard = tmp_path / "dashboard.html"
    rows = "".join(
        f"<tr><td>GENE{i}</td><td>guide{i}</td><td>{i / 100}</td><td>{i / 1000}</td></tr>"
        for i in range(40)
    )
    final_dashboard.write_text(
        '<h3 id="local">Guide Inference</h3><p>Local Analysis</p><table>'
        '<thead><tr><th>gene_name</th><th>guide_id</th><th>sceptre_log2_fc</th>'
        f'<th>sceptre_p_value</th></tr></thead><tbody>{rows}</tbody></table>',
        encoding="utf-8",
    )

    result = dashboard.final_inference_content(final_dashboard, limit=25)

    assert "Local Analysis: top guide–gene pairs" in result
    assert "Filter the 25 mirrored rows" in result
    assert "GENE24" in result
    assert "GENE25" not in result
    assert 'id="inference-0"' in result


def test_evaluation_skip_reason_and_artifacts_are_visible(tmp_path):
    evaluation = tmp_path / "evaluation_output"
    evaluation.mkdir()
    (evaluation / "controls_evaluation_skipped.txt").write_text(
        "Benchmark disabled because no validation set was supplied.\n", encoding="utf-8"
    )
    (evaluation / "local_analysis_sceptre.bedpe").write_text("chr8\t1\t2\n", encoding="utf-8")

    result = dashboard.evaluation_artifact_content(tmp_path)

    assert "Benchmark disabled" in result
    assert "local_analysis_sceptre.bedpe" in result
    assert "BEDPE" in result
