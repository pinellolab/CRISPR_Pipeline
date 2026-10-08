FROM condaforge/miniforge3:24.11.3-2

# Keep numba cache writable
ENV NUMBA_CACHE_DIR=/tmp

# Workdir and bring the env spec alongside
WORKDIR /app
COPY nextflow.yaml /app/

# (Optional) OS updates or packages; keep image slim
RUN apt-get update && \
    apt-get -y upgrade && \
    rm -rf /var/lib/apt/lists/*

# Create the env from your nextflow.yaml
# (use absolute path so it works regardless of WORKDIR)
RUN mamba env create -f /app/nextflow.yaml

# Fail the image build early if the scientific stack required by clone
# filtering and saturation QC cannot be imported together.
RUN /opt/conda/envs/nfenv/bin/python -c \
    "import anndata, matplotlib, mudata, numpy, pandas, scipy"

# Exercise the Scanpy sparse QC path during the build. This catches package
# resolutions that import successfully but fail only when Numba compiles QC.
RUN /opt/conda/envs/nfenv/bin/python -c \
    "import anndata as ad, numba, numpy as np, scanpy as sc; from scipy import sparse; assert tuple(map(int, numba.__version__.split('.')[:2])) < (0,67); a=ad.AnnData(sparse.csr_matrix(np.eye(501, dtype=np.int32))); sc.pp.calculate_qc_metrics(a, inplace=True, log1p=True, percent_top=(50,100,200,500)); assert 'pct_counts_in_top_500_genes' in a.obs"

# Put the env on PATH
# (assuming the env is named "nfenv" inside nextflow.yaml)
ENV PATH="/opt/conda/envs/nfenv/bin:${PATH}"

# # Copy any helper scripts/binaries (make sure build context includes ./bin)
# # If you don't have a bin/ in that folder, remove these two lines.
# COPY bin /usr/local/bin
# RUN chmod +x /usr/local/bin/* || true

# Default command
CMD ["/bin/bash", "--login"]
