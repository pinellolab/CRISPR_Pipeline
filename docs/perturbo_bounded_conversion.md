# PerTurbo result conversion memory

The adapter processes element and guide results separately, projects the columns
needed by the catalog, and writes bounded batches. The default batch size is
250,000 rows. Repeated identifiers use Arrow dictionary encoding.

BH correction still uses the complete applicable family: transcriptome-wide and
requested-pair outputs have separate corrections. It is never computed per
row group. Numeric family arrays remain in memory under a conservative
64-bytes-per-row planning budget, default 24 GiB. This budget is not a cap on
total process RSS: imports, mappings, I/O and conversion buffers also need RAM.
A family exceeding the configured budget fails before allocating its arrays.
As in the previous adapter, BH clips non-missing input probabilities to [0, 1]
(including infinities); raw output probabilities remain unchanged. All eight
available CRT diagnostic columns retain their source types and missing values.

The pipeline parameters are:

- `INFERENCE_PERTURBO_RESULT_BATCH_ROWS` (250000)
- `INFERENCE_PERTURBO_MAX_BH_WORKING_BYTES` (25769803776, 24 GiB)
- `INFERENCE_PERTURBO_COMPACT_RESULT_FLOATS` (false)

These do not increase the Nextflow task RAM allocation. Logs report requested
RAM, CPUs and container identity; available node/user RAM is not the task limit.
`MPLCONFIGDIR` points to a writable directory in task scratch.

## Precision

Linear p/q values and BH arithmetic remain float64. Current PerTurbo posterior
means, scales and probabilities are already float32, and the adapter generally
preserves those types. Optional compact storage rounds only allowlisted effect
and standard-error fields to float32, after inference/BH; original raw files
remain unchanged. It rejects finite overflow/underflow. Small rounding changes
are expected when the source effect fields were float64.

Catalogs currently retain raw linear p/q columns as well as derived log scores;
this change does not migrate the catalog to a log-only schema. The reusable
compaction helper also supports explicitly named derived log fields, but does
not change the existing catalog builders' defaults.

## Measured conversion memory

A CPU comparison on ml003 used saved essential-screen result prefixes, a
16-core affinity limit, and the existing runtime environment:

| Input rows | Old conversion peak RSS | New conversion peak RSS | Old/new seconds |
| --- | --- | --- | --- |
| 2,000,000 | 1.75 GiB | 516 MiB | 4.47 / 3.78 |
| 10,000,000 | not run | 754 MiB | — / 20.41 |

All 15 output columns matched exactly over two million rows, including p/q and
all eight CRT diagnostics. The reference was the existing converter on
`feat/standardize-inference-covariates` (18c09a0); both converters used the same
saved raw rows. The two-million-row peak RSS reduction was 71%.
This measures one table at a time; it does not include GPU inference. The
collaborator's 308-million-row input was unavailable, so that scale has not
been validated and these measurements are not a memory guarantee.

`--output-mudata` still explicitly materializes converted tables for embedding.
The production Nextflow process omits it by default; enabling it relinquishes
the bounded table-conversion memory behavior.

## Recovery

Raw fit outputs are created in task-local durable artifact scratch, not a
throwaway `/tmp` fit directory. A completed fit is published atomically with its
hash manifest before a sequential peer starts. Explicit parallel fitting still
works; successful peers are preserved if another fit fails.

Use `--conversion-only-artifact-dir PATH` with the original input path and fit
settings to regenerate tables after conversion failure. This verifies the raw
Parquets, producing metadata and compact saved identifier mappings. It does not
load or copy the full input MuData or launch GPU fitting. The input is stat-checked
for identity. These strict checks deliberately reject unrelated artifacts;
pre-existing outputs without the new metadata require reviewed manual adoption.
Derived requested-pair fallback tables belong to conversion scratch, leaving
raw published outputs immutable. Global and requested families remain distinct.

## Validation limits

The full Nextflow workflow was not launched here (Nextflow is unavailable in
this environment). Thirty-nine focused adapter/storage/merge tests passed at integration. They exercise output equivalence, empty families,
all-null early batches, recovery without data preparation, and interrupted
publication. Two existing `test_perturbo_result_parquet.py` failures reproduce
unchanged on dev d2f1ddd with the existing pandas runtime: categorical comparison
and mixed string/integer chromosome serialization. Those unrelated helpers were
not changed in this workstream. The standardized-covariate checks passed 22
tests, with one additional test blocked by a missing `seaborn` dependency.

## Suggested agent-guide addition

Future adapter changes should keep published native fit artifacts immutable,
validate recovery with real fit metadata present, and compare both global and
requested-family p/q values to the legacy converter. BH must not be performed
independently per batch. Storage dtype changes need separate validation from
scientific changes. These points would be useful additions to `AGENTS.md` when
this workflow is adopted.
