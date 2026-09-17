#!/usr/bin/env python3
"""Container-level smoke test for guide and post-clone HTO filters."""

import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).parents[1]


def make_mudata(path):
    cells = pd.Index([f"cell{i}" for i in range(8)])
    obs = pd.DataFrame({"batch": ["A"] * 6 + ["B"] * 2}, index=cells)
    gene = ad.AnnData(sparse.csr_matrix(np.ones((8, 2))), obs=obs.copy(), var=pd.DataFrame(index=["g1", "g2"]))
    guide = ad.AnnData(sparse.csr_matrix(np.ones((8, 3))), obs=obs.copy(), var=pd.DataFrame(index=["sg1", "sg2", "sg3"]))
    guide.layers["guide_assignment"] = sparse.csr_matrix(
        [[1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 1, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    hashing_obs = obs.copy()
    hashing_obs["hto_type_split"] = ["H1", "H1", "H2", "negative", "multiplets", "H2", "H3", "H3"]
    hashing = ad.AnnData(sparse.csr_matrix(np.ones((8, 3))), obs=hashing_obs, var=pd.DataFrame(index=["H1", "H2", "H3"]))
    md.MuData({"gene": gene, "guide": guide, "hashing": hashing}).write_h5mu(path)


def main():
    with tempfile.TemporaryDirectory(prefix="guide_hto_smoke_") as directory:
        tmp = Path(directory)
        source = tmp / "input.h5mu"
        make_mudata(source)
        guide_output = tmp / "guide_filtered.h5mu"
        subprocess.run(
            [
                "python", str(ROOT / "bin/filter_guide_assignment_qc.py"), str(source), str(guide_output),
                "--outdir", str(tmp / "guide_qc"), "--max-guides-per-cell", "2",
            ],
            check=True,
        )
        guide_filtered = md.read_h5mu(guide_output)
        assert guide_filtered.n_obs == 7 and "cell2" not in guide_filtered.obs_names

        # Emulate clone removal after guide QC. HTO support must be recomputed
        # only from this downstream population, not from the original cells.
        post_clone = guide_filtered[["cell0", "cell1", "cell3", "cell6", "cell7"]].copy()
        post_clone_path = tmp / "post_clone.h5mu"
        post_clone.write_h5mu(post_clone_path)
        hto_output = tmp / "hto_filtered.h5mu"
        subprocess.run(
            [
                "python", str(ROOT / "bin/filter_hto_post_clone.py"), str(post_clone_path), str(hto_output),
                "--outdir", str(tmp / "hto_qc"), "--min-positive-cells", "2", "--singlet-only", "true",
                "--filtered-hashing-output", str(tmp / "filtered.h5ad"),
                "--unfiltered-hashing-output", str(tmp / "unfiltered.h5ad"),
            ],
            check=True,
        )
        hto_filtered = md.read_h5mu(hto_output)
        assert hto_filtered.obs_names.tolist() == ["cell0", "cell1", "cell6", "cell7"]
        support = pd.read_csv(tmp / "hto_qc/hto_positive_cell_support.tsv", sep="\t")
        called = support.loc[support["called"], ["measurement_set", "hto_label"]]
        assert called.to_records(index=False).tolist() == [("A", "H1"), ("B", "H3")]
        print("guide/HTO QC smoke test passed")


if __name__ == "__main__":
    main()
