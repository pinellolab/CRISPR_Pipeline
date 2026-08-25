# Combined PerTurbo v2 JAX and cis-mask image

This build derives from Logan's immutable
`ghcr.io/pinellolab/perturbo:sha-39e62c1` image and adds the optional
`--gene-by-element-varm-key` and
`--gene-by-element-names-uns-key` local/cis pair mask.

The resulting single image supports both pipeline modes:

- local/cis inference passes the pair-mask arguments;
- global/trans inference omits them and retains the original Logan behavior.

Published image:

```text
ghcr.io/pinellolab/crispr_pipeline/perturbo:v2-cis-mask
```

The workflow validates the embedded Logan revision, the two mask CLI options,
the mask-specific regression tests, and the unmasked output path before push.
The source revisions are Logan/JAX `39e62c18a64e73b3da3b2c5dba12fda1e7d62e52`
plus cis-mask patch `9bfcc052231d30e97f21d9ecac411aab52a23938`.
