import pathlib
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sparse = pytest.importorskip("scipy.sparse")

import anndata_concat


def _write_input(tmp_path, batch_num, x, layers=None, folder="counts_unfiltered"):
    """Build the <dir>_ks_<suffix>/<folder>/adata.h5ad structure
    anndata_concat expects, with a directory name extract_batch_num can
    parse (regex: `(.+)_ks_`)."""
    n_obs, n_var = x.shape
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_var)]),
    )
    for name, values in (layers or {}).items():
        adata.layers[name] = values

    file_dir = tmp_path / f"{batch_num}_ks_out"
    counts_dir = file_dir / folder
    counts_dir.mkdir(parents=True)
    adata.write_h5ad(counts_dir / "adata.h5ad")
    return str(file_dir)


def _write_covariates(tmp_path, batch_num):
    covariate_path = tmp_path / "covariates.csv"
    pd.DataFrame({"batch": [str(batch_num)], "barcode_key": [str(batch_num)]}).to_csv(
        covariate_path, index=False
    )
    return str(covariate_path)


def _run_main(monkeypatch, args):
    monkeypatch.setattr(sys, "argv", ["anndata_concat.py"] + args)
    anndata_concat.main()


def test_standard_workflow_output_is_sparse_narrow_dtype(tmp_path, monkeypatch):
    """A plain dense float64 X (no nac layers) should end up sparse and
    narrowed to uint16 once every value is a whole, non-negative number."""
    x = np.array([[0.0, 3.0], [7.0, 0.0]], dtype=np.float64)
    file_dir = _write_input(tmp_path, batch_num=0, x=x)
    covariate_path = _write_covariates(tmp_path, batch_num=0)
    output_path = tmp_path / "combined.h5ad"

    _run_main(
        monkeypatch,
        [
            file_dir,
            covariate_path,
            "--output", str(output_path),
            "--temp_dir", str(tmp_path / "temp"),
        ],
    )

    result = ad.read_h5ad(output_path)
    assert sparse.issparse(result.X)
    assert result.X.dtype == np.uint16
    np.testing.assert_array_equal(result.X.toarray(), x)


def test_nac_workflow_stays_float_when_fractional_and_mm_disabled(tmp_path, monkeypatch):
    """cc-perturb-seq's nac (nascent/mature/ambiguous) workflow can produce
    fractional counts from EM-distributed ambiguous reads. With mm=false
    (no rounding), the combined .X must stay float, not get silently
    truncated to an integer dtype."""
    mature = np.array([[1.0, 2.0]])
    nascent = np.array([[0.0, 1.0]])
    ambiguous = np.array([[0.5, 0.25]])  # fractional -- the case that matters
    x_placeholder = mature + nascent + ambiguous  # what kb-python would put in X

    file_dir = _write_input(
        tmp_path,
        batch_num=0,
        x=x_placeholder,
        layers={"mature": mature, "nascent": nascent, "ambiguous": ambiguous},
    )
    covariate_path = _write_covariates(tmp_path, batch_num=0)
    output_path = tmp_path / "combined.h5ad"

    _run_main(
        monkeypatch,
        [
            file_dir,
            covariate_path,
            "--output", str(output_path),
            "--temp_dir", str(tmp_path / "temp"),
            "--mm", "false",
        ],
    )

    result = ad.read_h5ad(output_path)
    assert sparse.issparse(result.X)
    assert result.X.dtype == np.float32, "fractional nac counts must not be truncated to int"
    expected = mature + nascent + ambiguous
    np.testing.assert_allclose(result.X.toarray(), expected)
    # The individual layers are never rounded and must also stay float.
    assert result.layers["ambiguous"].dtype == np.float32
    np.testing.assert_allclose(result.layers["ambiguous"].toarray(), ambiguous)


def test_nac_workflow_narrows_to_int_when_mm_rounds_to_whole_numbers(tmp_path, monkeypatch):
    """With mm=true, .X gets rounded to whole numbers -- it should then be
    stored as a narrow int dtype, even though the pre-rounding sum was
    fractional. The un-rounded layers stay float regardless."""
    mature = np.array([[1.0, 2.0]])
    nascent = np.array([[0.0, 1.0]])
    ambiguous = np.array([[0.5, 0.25]])
    x_placeholder = mature + nascent + ambiguous

    file_dir = _write_input(
        tmp_path,
        batch_num=0,
        x=x_placeholder,
        layers={"mature": mature, "nascent": nascent, "ambiguous": ambiguous},
    )
    covariate_path = _write_covariates(tmp_path, batch_num=0)
    output_path = tmp_path / "combined.h5ad"

    _run_main(
        monkeypatch,
        [
            file_dir,
            covariate_path,
            "--output", str(output_path),
            "--temp_dir", str(tmp_path / "temp"),
            "--mm", "true",
        ],
    )

    result = ad.read_h5ad(output_path)
    assert sparse.issparse(result.X)
    assert result.X.dtype == np.uint16
    expected_rounded = np.round(mature + nascent + ambiguous)
    np.testing.assert_array_equal(result.X.toarray(), expected_rounded)
    # Layers are untouched by rounding and stay float since they're still
    # genuinely fractional.
    assert result.layers["ambiguous"].dtype == np.float32
