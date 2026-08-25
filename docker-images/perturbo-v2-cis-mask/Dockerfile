ARG LOGAN_PERTURBO_IMAGE=ghcr.io/pinellolab/perturbo:sha-39e62c1
FROM ${LOGAN_PERTURBO_IMAGE}

LABEL org.opencontainers.image.title="PerTurbo v2 JAX with cis mask" \
      org.opencontainers.image.description="Logan PerTurbo v2 JAX plus optional local gene-by-element pair masking" \
      org.opencontainers.image.source="https://github.com/pinellolab/CRISPR_Pipeline" \
      io.pinellolab.perturbo.logan-revision="39e62c18a64e73b3da3b2c5dba12fda1e7d62e52" \
      io.pinellolab.perturbo.cis-mask-revision="9bfcc052231d30e97f21d9ecac411aab52a23938"

COPY perturbo-cis-mask.patch /tmp/perturbo-cis-mask.patch

# The Logan image keeps its source under /app/src and installs a non-editable
# copy under /usr/local. Patch both representations so the CLI and imported
# package always use the same combined implementation.
RUN set -eu; \
    apt-get update; \
    apt-get install -y --no-install-recommends patch; \
    cd /app; \
    patch --batch --forward -p1 < /tmp/perturbo-cis-mask.patch; \
    package_dir="$(python -c 'from pathlib import Path; import perturbo; print(Path(perturbo.__file__).parent)')"; \
    cp /app/src/perturbo/cli.py "$package_dir/cli.py"; \
    cp /app/src/perturbo/core.py "$package_dir/core.py"; \
    cp /app/src/perturbo/results.py "$package_dir/results.py"; \
    find "$package_dir" -type f -name '*.pyc' -delete; \
    rm -f /tmp/perturbo-cis-mask.patch; \
    apt-get purge -y --auto-remove patch; \
    rm -rf /var/lib/apt/lists/*

RUN set -eu; \
    command -v perturbo; \
    perturbo --help | grep -F -- '--gene-by-element-varm-key'; \
    perturbo --help | grep -F -- '--gene-by-element-names-uns-key'; \
    python -c 'import perturbo, pyarrow'

CMD ["perturbo", "--help"]
