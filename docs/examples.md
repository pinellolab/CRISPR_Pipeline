# Example configs and datasets

These examples target **CRISPR Pipeline v1.1.0**. Install the bundled dashboard
plugin, pin Nextflow 26.04.6, and use the v1 parser as described in
[Slurm setup](run-recipes.md#slurm-configuration).

## Configuration downloads

Pass one assay overlay with `-c` on top of the pipeline defaults. Supply
the samplesheet, output directory, and site settings from your run.

| Config | Use |
|---|---|
| [nextflow_cc.config](downloadable-examples/nextflow_cc.config) | CC-Perturb-seq, including barcode replacement. Check the barcode asset location and reference annotation. |
| [nextflow_tapseq.config](downloadable-examples/nextflow_tapseq.config) | Targeted TAP-seq. Supply the panel GTF via `--REFERENCE_gtf_local_path`. |
| [nextflow_yang.dev.config](downloadable-examples/nextflow_yang.dev.config) | Compact Yang assay and Slurm overlay. Defaults to a full run. |
| [nextflow_slurm_all_options.config](downloadable-examples/nextflow_slurm_all_options.config) | Commented engine defaults and optional settings grouped by topic. Uncomment only the intended overrides. |

## Study inputs

Each ZIP contains the v1.1.0-compatible assay config, samplesheet, measurement-set
metadata, guide metadata, and modality SeqSpecs. The Gary Hon example also
includes HTO metadata and a FASTQ accession manifest. Get the reads and cell-barcode
allowlist from the sources in the dataset README. The file layout in each
samplesheet preserves the corresponding run's modality and lane pairing.

| Dataset | Measurements | Downloads |
|---|---:|---|
| Gasperini at-scale K562 CRISPRi | 32 sets · 388 rows | [ZIP](downloadable-examples/gasperini_at_scale.zip) · [Config](downloadable-examples/gasperini_at_scale/nextflow.config) · [Samplesheet](downloadable-examples/gasperini_at_scale/samplesheet.csv) · [Input notes](downloadable-examples/gasperini_at_scale/README.md) |
| Replogle K562 day-6 essential genes | 48 GEM groups · 576 rows | [ZIP](downloadable-examples/replogle_k562_essential.zip) · [Config](downloadable-examples/replogle_k562_essential/nextflow.config) · [Samplesheet](downloadable-examples/replogle_k562_essential/samplesheet.csv) · [Input notes](downloadable-examples/replogle_k562_essential/README.md) |
| Gary Hon WTC11 TF benchmark | 3 sets · 12 rows | [ZIP](downloadable-examples/gary_hon_tf_benchmark.zip) · [Config](downloadable-examples/gary_hon_tf_benchmark/nextflow.config) · [Samplesheet](downloadable-examples/gary_hon_tf_benchmark/samplesheet.csv) · [Input notes](downloadable-examples/gary_hon_tf_benchmark/README.md) |

The [examples directory on GitHub](https://github.com/pinellolab/CRISPR_Pipeline/tree/main/examples)
also contains these files. Use the [tagged pipeline checkout](https://github.com/pinellolab/CRISPR_Pipeline/tree/v1.1.0)
to execute them.

### Resolve the samplesheet

Each ZIP includes `prepare_samplesheet.py`. Run it against the supplied
samplesheet with the local directories containing the reads and allowlist.
The script resolves the included metadata and SeqSpecs inside the example
directory, verifies that every input exists, and writes a new samplesheet.

```bash
python3 /path/to/example/prepare_samplesheet.py /path/to/example/samplesheet.csv \
  --fastq-dir /shared/data/reads --barcode-dir /shared/data/allowlists \
  --output /shared/run/samplesheet.csv

nextflow -c /path/to/example/nextflow.config -c /path/to/site.config \
  run main.nf -profile slurm \
  --input /shared/run/samplesheet.csv --outdir /shared/run/results \
  --ANALYSIS_mode checkpoint -work-dir /shared/run/work
```

Run Nextflow from the release checkout. Set your site's Slurm queue, container
runtime and shared cache using the [Slurm recipe](run-recipes.md#slurm-configuration).
Use `--ANALYSIS_mode full` for inference in the same invocation. The checkpoint
is the recommended first run for input and cohort review.

## Input formats

For smaller format references, see the [minimal RNA/guide samplesheet](downloadable-examples/samplesheet.csv),
[hashing samplesheet](downloadable-examples/hashing_samplesheet.csv),
[RNA SeqSpec](downloadable-examples/rna_seqspec.yml),
[guide SeqSpec](downloadable-examples/guide_seqspec.yml),
[HTO SeqSpec](downloadable-examples/hash_seqspec.yml),
[guide metadata](downloadable-examples/guide_metadata.tsv), and
[HTO metadata](downloadable-examples/hash_metadata.tsv).

## Completed benchmark output

The [Gary Hon full-depth benchmark on Synapse](https://www.synapse.org/Synapse:syn78349990)
contains the completed inference MuData, Parquet tables, QC archive, dashboard,
and run provenance. See the [results gallery](output-gallery.md) for plots
and the [output guide](output.md) for file definitions.
