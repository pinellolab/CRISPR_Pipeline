# CRISPR Pipeline examples

These inputs and overlays target **v1.1.0**. See the
[website's examples guide](https://pinellolab.github.io/CRISPR_Pipeline/examples/)
for downloads and launch instructions.

- [Gasperini at-scale K562 CRISPRi](datasets/gasperini_at_scale/)
- [Replogle K562 day-6 essential genes](datasets/replogle_k562_essential/)
- [Gary Hon WTC11 TF benchmark](datasets/gary_hon_tf_benchmark/)

Supply FASTQs and the cell-barcode allowlists separately. From the repository
root, run `python3 examples/prepare_samplesheet.py --help` to localize a
samplesheet. Each assay config inherits defaults and containers from the
pipeline root. Site profiles and resource settings are documented on the site.
