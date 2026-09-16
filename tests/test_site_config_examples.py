"""The example site configs must be deltas from nextflow.config, and nothing else.

nextflow.config is always loaded, even when an example is passed with `-c`, so a
line in an example that merely restates a pipeline default does nothing today
and silently disagrees with nextflow.config the moment the default moves --
which is what happened to INFERENCE_PERTURBO_MAX_CHUNK_CELLS and
containers.perturbo before the header comments went in. These tests enforce that
invariant mechanically instead of by comment, and pin the TAP-seq example's
panel-defining parameters.
"""

import pathlib
import re

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BASE_CONFIG = REPO_ROOT / "nextflow.config"
# Only the example this change adds. The CC-Perturb-seq config belongs to that
# assay: it carries a known restatement of a nextflow.config default, and
# asserting on it from here would fail our suite for an edit that is not ours.
SITE_CONFIGS = ("nextflow_tapseq.config",)

# Deltas a site config is allowed to state even when they match the pipeline
# default: the run's own inputs, and the ceilings every site must look at.
EXEMPT = frozenset({"input", "outdir", "max_cpus", "max_memory"})


def _strip_comment(line):
    """Drop a trailing `//` comment without cutting inside a quoted string."""
    out = []
    quote = None
    i = 0
    while i < len(line):
        ch = line[i]
        if quote is not None:
            out.append(ch)
            if ch == "\\" and i + 1 < len(line):
                out.append(line[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "/" and line[i + 1:i + 2] == "/":
            break
        else:
            out.append(ch)
        i += 1
    return "".join(out)


def _mask_strings(line):
    """Blank out string contents so brace counting cannot trip over `${...}`."""
    out = []
    quote = None
    for ch in line:
        if quote is not None:
            out.append("_" if ch != quote else ch)
            if ch == quote:
                quote = None
        else:
            out.append(ch)
            if ch in "'\"":
                quote = ch
    return "".join(out)


def _normalize(value):
    """Compare 1e-06 with 0.000001, and '128.GB' with 128.GB, as equal."""
    value = value.strip().rstrip(";")
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
        value = value[1:-1]
    if value in ("true", "false"):
        return value == "true"
    if value == "null":
        return None
    try:
        return float(value)
    except ValueError:
        return value


_ASSIGN = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\S.*)$")


def _params(path):
    """Top-level assignments in every `params { }` block, last one winning.

    nextflow.config carries three params blocks (the Pipeline Configurator
    appends two), and Nextflow applies them in order, so the effective default
    is the last assignment -- QC_barcode_filter's real default is the injected
    'knee', not the value in the first block.
    """
    values = {}
    depth = 0
    params_depth = None
    for raw in path.read_text().splitlines():
        code = _strip_comment(raw).strip()
        masked = _mask_strings(code)
        if params_depth is not None and depth == params_depth and not masked.endswith("{"):
            match = _ASSIGN.match(code)
            if match:
                values[match.group(1)] = _normalize(match.group(2))
        if params_depth is None and re.match(r"^params\s*\{$", masked):
            params_depth = depth + 1
        depth += masked.count("{") - masked.count("}")
        if params_depth is not None and depth < params_depth:
            params_depth = None
    return values


def test_the_parser_sees_the_defaults_it_is_asked_about():
    """Guard the guard: a parser that silently found nothing would pass anything."""
    defaults = _params(BASE_CONFIG)
    assert defaults["TAPSEQ_QC_MODE"] is False
    assert defaults["QC_min_genes_per_cell"] == 500
    assert defaults["QC_min_cells_per_gene"] == 0.05
    assert defaults["QC_barcode_filter"] == "knee"  # from the last params block
    assert defaults["GUIDE_ASSIGNMENT_capture_method"] == "direct-capture"
    assert defaults["spacer_tag"] == "TAGCTCTTAAAC"
    assert defaults["reverse_complement_guides"] is True
    assert defaults["REFERENCE_restrict_genes_to_gtf"] is False
    assert defaults["INFERENCE_control_group"] == "auto"
    assert defaults["INFERENCE_PERTURBO_CRT_POOL"] == "from-moi"
    assert defaults["INFERENCE_SCEPTRE_GENE_CHUNK_SIZE"] == 1000
    assert "containers" not in defaults  # nested block, not an assignment


def _restated_defaults(name):
    defaults = _params(BASE_CONFIG)
    return {
        key: value
        for key, value in _params(REPO_ROOT / name).items()
        if key not in EXEMPT and key in defaults and defaults[key] == value
    }


@pytest.mark.parametrize("name", SITE_CONFIGS)
def test_a_site_config_states_only_deltas(name):
    restated = _restated_defaults(name)
    assert not restated, (
        f"{name} restates nextflow.config defaults, which will drift: {restated}. "
        "Delete the lines; nextflow.config is loaded either way."
    )


@pytest.mark.parametrize("name", SITE_CONFIGS)
def test_a_site_config_does_not_repin_the_containers(name):
    """The pins live in one place. An example that copies them goes stale, and
    containers.perturbo is the case that actually did."""
    text = (REPO_ROOT / name).read_text()

    assert not re.search(r"^\s*containers\s*\{", text, re.M), name
    assert "ghcr.io" not in text, name
    assert not re.search(r"^\s*perturbo\s*=", text, re.M), name


@pytest.mark.parametrize("name", SITE_CONFIGS)
def test_a_site_config_says_why_it_only_carries_deltas(name):
    text = (REPO_ROOT / name).read_text()
    assert "always loaded automatically" in text, name
    assert "DELTA" in text, name
    assert "includeConfig 'conf/provenance.config'" in text, name


@pytest.mark.parametrize("name", SITE_CONFIGS)
def test_a_site_config_is_discoverable_from_the_readme(name):
    assert name in (REPO_ROOT / "README.md").read_text(), name


def test_the_tapseq_example_defines_the_panel():
    """The four settings that make a targeted screen testable at all."""
    params = _params(REPO_ROOT / "nextflow_tapseq.config")

    assert params["TAPSEQ_QC_MODE"] is True
    assert params["REFERENCE_restrict_genes_to_gtf"] is True
    # Whole-transcriptome floors: a cell expressing 60 of 68 panel genes must pass.
    assert params["QC_min_genes_per_cell"] < 500
    # Redundant once the panel restriction has dropped the untestable genes.
    assert params["QC_min_cells_per_gene"] < 0.05


def test_the_tapseq_example_leaves_the_panel_gtf_to_the_user():
    """REFERENCE_gtf_local_path must be the screen's panel, so it cannot be shipped.

    Hard-coding a path here would be worse than leaving it out: the pipeline
    falls back to the full GENCODE download when the local path is missing, and
    a transcriptome-wide GTF makes REFERENCE_restrict_genes_to_gtf a no-op.
    """
    text = (REPO_ROOT / "nextflow_tapseq.config").read_text()

    assert "REFERENCE_gtf_local_path" not in _params(REPO_ROOT / "nextflow_tapseq.config")
    assert re.search(r"^\s*//\s*REFERENCE_gtf_local_path\s*=", text, re.M), (
        "the required user value must still be visible as a commented TO-DO"
    )
    assert "TO-DO" in text


def test_the_tapseq_example_leaves_the_control_group_shared():
    """Neither the shared setting nor a per-method override is set.

    A per-method parameter overrides the shared setting for one method only, so
    an example that set INFERENCE_PERTURBO_CRT_POOL would silently make SCEPTRE
    and PerTurbo answer different questions on every screen that copied it.
    """
    params = _params(REPO_ROOT / "nextflow_tapseq.config")

    for key in ("INFERENCE_control_group", "INFERENCE_PERTURBO_CRT_POOL",
                "INFERENCE_SCEPTRE_control_group", "Multiplicity_of_infection"):
        assert key not in params, f"{key} must stay at nextflow.config's default"


def test_the_tapseq_example_records_the_single_lane_lesson():
    """`auto` is wrong when the non-targeting cells are one batch, and that is
    invisible from the MOI. The example has to say so, or it teaches the bug."""
    text = (REPO_ROOT / "nextflow_tapseq.config").read_text()

    assert "INFERENCE_control_group" in text
    assert "complement" in text
    assert re.search(r"\bcrosstab\b", text, re.I), "tell the reader how to check their own screen"
    assert re.search(r"lane|batch", text, re.I)


def test_the_tapseq_example_marks_the_library_chemistry_as_dataset_specific():
    """These four worked for chr8; copying them unverified is the failure mode."""
    params = _params(REPO_ROOT / "nextflow_tapseq.config")
    text = (REPO_ROOT / "nextflow_tapseq.config").read_text()

    assert params["GUIDE_ASSIGNMENT_capture_method"] == "crop-seq"
    assert params["spacer_tag"] == ""
    assert params["reverse_complement_guides"] is False
    assert params["QC_barcode_filter"] == "knee2"
    assert "seqspec" in text.lower(), "point the reader at their own seqspec check"
