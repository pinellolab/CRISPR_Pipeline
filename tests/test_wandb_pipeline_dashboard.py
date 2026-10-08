import argparse
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "render_wandb_pipeline_dashboard.py"
SPEC = importlib.util.spec_from_file_location("wandb_pipeline_dashboard", SCRIPT)
dashboard = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(dashboard)


def test_special_qc_processes_are_routed_to_expected_categories():
    assert dashboard.family_for("workflow:skipGTFDownload") == "input"
    assert dashboard.family_for("workflow:sequencing_saturation") == "evaluation"
    assert dashboard.family_for("workflow:remove_clonal_cells") == "postconcat_qc"
    assert dashboard.image_family(Path("measurement_set_qc/rna_qc_filter_flow_B1.png")) == "preprocessing"


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
    assert 'for="family-tab-seqspec"' in result
    assert 'id="family-tab-seqspec"' in result
    assert 'id="family-inference"' in result
    assert 'id="family-postconcat_qc"' in result
    assert "FAILED" in result
    assert "42" in result
    assert "data:image/png;base64," in result
    assert "No credentials, FASTQs or unbounded task logs embedded" in result
    assert 'type="radio" name="pipeline-family"' in result
    assert "familyFromHash" not in result
    assert '#family-tab-seqspec:checked~.family-panels #family-seqspec{display:block}' in result


def test_processes_are_mapped_to_expected_families():
    assert dashboard.family_for("pipeline:seqspeccheck") == "seqspec"
    assert dashboard.family_for("pipeline:mapping_rna_pipeline:mappingrna") == "mapping"
    assert dashboard.family_for("pipeline:guide_assignment_sceptre") == "guide_assignment"
    assert dashboard.family_for("pipeline:inference_pipeline:inference_perturbo") == "inference"
    assert dashboard.family_for("pipeline:additional_qc") == "evaluation"
    assert dashboard.family_for("pipeline:dashboard") == "final"


def test_live_tasks_are_added_from_nextflow_log_without_false_failure_panel(tmp_path):
    trace_rows = [{"task_id": "1", "process": "pipeline:mappingGuide", "status": "COMPLETED"}]
    log = tmp_path / "nextflow.log"
    log.write_text(
        "DEBUG Unable to get file attributes -- Cause: missing\n"
        "TaskHandler[id: 2; name: NFCORE_CRISPR:CRISPR_PIPELINE:mapping_rna_pipeline:mappingscRNA (1); "
        "status: RUNNING; exit: -; error: -; workDir: /work/aa/bb]\n",
        encoding="utf-8",
    )
    merged = dashboard.merge_live_tasks(trace_rows, log)
    assert len(merged) == 2
    assert merged[-1]["status"] == "RUNNING"
    assert dashboard.family_for(merged[-1]["process"]) == "mapping"
    assert dashboard.failure_content(merged, log, 30) == ""


def test_newly_submitted_task_is_visible_before_taskhandler_update(tmp_path):
    log = tmp_path / "nextflow.log"
    log.write_text(
        "Submitted process > NFCORE_CRISPR:CRISPR_PIPELINE:mapping_rna_pipeline:mappingscRNA (4)\n",
        encoding="utf-8",
    )
    merged = dashboard.merge_live_tasks([], log)
    assert len(merged) == 1
    assert merged[0]["status"] == "SUBMITTED"
    assert dashboard.family_for(merged[0]["process"]) == "mapping"


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


def test_mapping_cards_group_modalities_by_measurement_set():
    data = {
        "observed_metrics": {
            "mapping_json": [
                {
                    "measurement_set": "SET_A",
                    "modality": "RNA",
                    "metrics": {
                        "n_processed": 1000,
                        "p_pseudoaligned": 88.5,
                        "p_unique": 80.0,
                    },
                },
                {
                    "measurement_set": "SET_A",
                    "modality": "guide",
                    "metrics": {
                        "n_processed": 500,
                        "p_pseudoaligned": 75.0,
                        "percentageReadsOnOnlist": 92.0,
                    },
                },
            ]
        }
    }

    result = dashboard.mapping_measurement_set_cards(data)

    assert "Mapping by measurement set" in result
    assert result.count('<details class="measurement-card">') == 1
    assert "SET_A" in result
    assert "RNA" in result and "guide" in result
    assert "1,500" in result
    assert "88.5%" in result


def test_preprocessing_cards_keep_each_sets_plots_inside_its_card():
    first_knee = Path("knee_plot_scRNA_SET_A.png")
    first_flow = Path("rna_qc_filter_flow_SET_A.png")
    other_plot = Path("knee_plot_scRNA_SET_B.png")
    data = {
        "observed_metrics": {
            "measurement_set_rna_qc": {
                "rows": [
                    {
                        "measurement_set": "SET_A",
                        "input_barcodes": 1000,
                        "post_knee_cells": 300,
                        "fixed_min_counts": 500,
                        "post_min_counts_cells": 250,
                        "mad_total_counts_n": 5,
                        "post_mad_cells": 225,
                        "scrublet_enabled": True,
                        "removed_by_scrublet": 25,
                        "retained_cells": 200,
                        "knee_umi_threshold": 42,
                        "knee_rank": 300,
                    }
                ]
            }
        }
    }

    result, used = dashboard.preprocessing_measurement_set_cards(
        data,
        [(first_knee, "aW1hZ2U="), (first_flow, "Zmxvdw=="), (other_plot, "b3RoZXI=")],
    )

    assert "Preprocessing by measurement set" in result
    assert "automatic knee → UMI → MAD → doublet policy" in result
    assert "After Scrublet" in result
    assert "20.00%" in result
    assert "knee plot scRNA SET A" in result
    assert "rna qc filter flow SET A" in result
    assert "SET_B" not in result
    assert used == {first_knee, first_flow}


def test_absolute_output_name_does_not_force_every_image_into_postconcat():
    root = Path("/results/results_embedding_qc_20260921")
    assert dashboard.image_family(root / "measurement_set_qc/figures/knee_plot_scRNA_SET_A.png") == "preprocessing"
    assert dashboard.image_family(root / "guide_assignment_qc/guide_assignment_filter_steps_SET_A.png") == "guide_assignment"
    assert dashboard.image_family(root / "pipeline_dashboard/additional_qc/intended_target_volcano.png") == "evaluation"
    assert dashboard.image_family(root / "postconcat_embedding_qc/after_clone/embedding_qc/pca_qc_panel.png") == "postconcat_qc"


def test_figure_cards_group_stage_plots_by_measurement_set():
    first = Path("guide_assignment_filter_steps_SET_A.png")
    second = Path("guide_assignment_filter_steps_SET_B.png")

    result, used = dashboard.figure_measurement_set_cards(
        [(first, "Zmlyc3Q="), (second, "c2Vjb25k")],
        ["SET_A", "SET_B"],
        "Guide assignment by measurement set",
        "Open a card.",
        "guide calling",
    )

    assert "Guide assignment by measurement set" in result
    assert result.count('<details class="measurement-card figure-card">') == 2
    assert "SET_A" in result and "SET_B" in result
    assert used == {first, second}


def test_completed_dashboard_opens_on_measurement_set_hierarchy(tmp_path):
    trace = tmp_path / "trace.tsv"
    trace.write_text("process\tname\tstatus\nworkflow:PreprocessAnnData\tpreprocess\tCOMPLETED\n", encoding="utf-8")
    guide = tmp_path / "guides.json"
    guide.write_text("{}", encoding="utf-8")
    seqspec = tmp_path / "seqspec.csv"
    seqspec.write_text("Sample,IsWinner\n", encoding="utf-8")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        '{"observed_metrics":{"measurement_set_rna_qc":{"rows":['
        '{"measurement_set":"SET_A","input_barcodes":100,"retained_cells":80}]}}}',
        encoding="utf-8",
    )

    result = dashboard.render(
        argparse.Namespace(
            trace=trace,
            run_id="completed-run",
            run_name="completed run",
            status="completed",
            guide_report=guide,
            seqspec_table=seqspec,
            seqspec_image=None,
            qc_metrics_json=metrics,
            artifact_dir=None,
            final_dashboard_html=None,
            nextflow_log=None,
            max_image_bytes=1_000_000,
            tail_lines=10,
        )
    )

    assert 'for="family-tab-input"' in result
    assert 'for="family-tab-inference"' in result
    assert 'id="family-tab-preprocessing" checked' in result
    assert 'id="family-preprocessing" class="family-panel"' in result
    assert 'onclick="selectFamily' not in result
    assert "selectFamily" not in result
