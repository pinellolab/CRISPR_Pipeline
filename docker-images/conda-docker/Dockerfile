FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV NUMBA_CACHE_DIR=/tmp
ENV PATH="/app/.venv/bin:${PATH}"

WORKDIR /app
COPY requirements.txt /app/

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    curl \
    git \
    default-jre-headless && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow as a standalone CLI.
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/nextflow && \
    chmod +x /usr/local/bin/nextflow

# Install Python dependencies with uv into an isolated venv.
RUN pip install --no-cache-dir uv && \
    uv venv /app/.venv && \
    uv pip install --python /app/.venv/bin/python \
    -r /app/requirements.txt

CMD ["/bin/bash"]
