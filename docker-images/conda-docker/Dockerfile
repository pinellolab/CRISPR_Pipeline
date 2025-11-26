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

# Put the env on PATH
# (assuming the env is named "nfenv" inside nextflow.yaml)
ENV PATH="/opt/conda/envs/nfenv/bin:${PATH}"

RUN mamba run -n nfenv mamba upgrade -y scanpy
RUN mamba run -n nfenv pip install --no-cache-dir \
  "git+https://github.com/IGVF-DACC/seqspec.git@v25-09-23"

# # Copy any helper scripts/binaries (make sure build context includes ./bin)
# # If you don't have a bin/ in that folder, remove these two lines.
# COPY bin /usr/local/bin
# RUN chmod +x /usr/local/bin/* || true

# Symlink kallisto & bustools provided by kb-python (path may vary with Python version)
# If Python version differs, adjust the 3.12 segment accordingly.
RUN ln -sf /opt/conda/envs/nfenv/lib/python3.12/site-packages/kb_python/bins/linux/*/kallisto /opt/conda/envs/nfenv/bin/ && \
    ln -sf /opt/conda/envs/nfenv/lib/python3.12/site-packages/kb_python/bins/linux/*/bustools /opt/conda/envs/nfenv/bin/

# Default command
CMD ["/bin/bash", "--login"]
