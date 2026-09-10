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
    assert "No credentials, FASTQs or task logs embedded" in result


def test_processes_are_mapped_to_expected_families():
    assert dashboard.family_for("pipeline:seqspeccheck") == "seqspec"
    assert dashboard.family_for("pipeline:mapping_rna_pipeline:mappingrna") == "mapping"
    assert dashboard.family_for("pipeline:guide_assignment_sceptre") == "guide_assignment"
    assert dashboard.family_for("pipeline:inference_pipeline:inference_perturbo") == "inference"
    assert dashboard.family_for("pipeline:additional_qc") == "evaluation"
    assert dashboard.family_for("pipeline:dashboard") == "final"
