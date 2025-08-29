
# IGVF single-cell Perturb-seq Pipeline Data Downloader

This toolkit provides a two-step workflow to generate per-sample metadata and download the corresponding files from the [IGVF data portal](https://data.igvf.org) for Perturb-seq pipelines.

---

## Requirements

- Python 3.7+
- IGVF API credentials (keypair JSON)
- (Optional) Google Cloud credentials (`gcloud auth application-default login` OR set environment variable for `GOOGLE_APPLICATION_CREDENTIALS`="/path/to/your-service-account.json")
- Install dependencies:

```bash
pip install -r requirements.txt
```

---

## API Key Format

Create a file `igvf_key.json` with the following structure:

```json
{
  "key": "IGVF_API_KEY",
  "secret": "IGVF_SECRET_KEY"
}
```

OR set environment variables: `IGVF_API_KEY` and `IGVF_SECRET_KEY`.

---

## Step 1: Generate Per-Sample Metadata TSV

### If `seqspec` files are available in the IGVF data portal

```bash
python3 generate_per_sample.py \
  --keypair igvf_key.json \
  --accession IGVFDS9445RJOU \
  --output test_fetch.tsv
```

This will generate a `test_fetch.tsv` file that contains a per-sample metadata sheet with file accessions and experimental details.

### If `seqspec` files are **missing**, provide fallback YAMLs

```bash
python3 generate_per_sample.py \
  --keypair igvf_key.json \
  --accession IGVFDS7340YDHF \
  --output test_fetch.tsv \
  --hash_seqspec hash_seq_spec.yaml \
  --rna_seqspec rna_seq_spec.yaml \
  --sgrna_seqspec sgrna_seq_spec.yaml
```

Fallbacks are required per modality:
- `--rna_seqspec`: scRNA-seq
- `--sgrna_seqspec`: sgRNA
- `--hash_seqspec`: cell hashing

---

## Step 2: Download and Verify Files

Use the `.tsv` output from Step 1 as input here.

### Download All Files (FASTQ + Reference)

```bash
python3 download_igvf.py \
  --sample test_fetch.tsv \
  --gunzip \
  --keypair igvf_key.json
```

This will:
- Download FASTQ, YAML, and TSV files from the IGVF data portal
- Verify MD5 checksums
- Gunzip files and store them locally
- Generate `updated_paths_test_fetch.tsv` with local file paths

---

## Upload to Google Cloud Storage (Optional)

### Upload files to GCS while downloading:

```bash
python3 download_igvf.py \
  --sample test_fetch.tsv \
  --keypair igvf_key.json \
  --gcs-upload \
  --gcs-bucket google-cloud-bucket \
  --gcs-prefix pipeline_runs/IGVFDS7340YDHF_2025_08_08/ \
  --output-dir ./downloads
```

This will:
- Upload each file to the specified bucket and prefix
- Update the sample sheet to use `gs://` paths
- Upload the updated sample sheet to the bucket

### Prevent overwrite of existing GCS files (default)

If you want to **force upload** to the same path on google cloud, add:

```bash
--gcs-force-upload
```

---

## Dry Run

Preview download/upload actions without performing them:

```bash
--dry-run
```

---

## Optional Arguments

| Argument              | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `--file-types fastq`  | Download only FASTQ files                                                   |
| `--file-types other`  | Download only YAML/TSV configuration files                                  |
| `--file-types all`    | *(default)* Download everything                                             |
| `--output-dir`        | Output directory for locally downloaded files (default: `./downloads`)             |
| `--gunzip`        | Gunzip downloaded files (default: keep compressed).             |
| `--gcs-upload`        | Enable GCS upload                                                           |
| `--gcs-bucket`        | Name of the GCS bucket                                                      |
| `--gcs-prefix`        | Subfolder path inside the GCS bucket (e.g. `pipeline_runs/IGVFDS...`)      |
| `--gcs-force-upload`  | Force overwrite on GCS if file already exists                               |
| `--dry-run`           | Print planned actions, but skip downloads and uploads                       |

---

## Output Files

- `test_fetch.tsv`: Metadata TSV from `generate_per_sample.py`
- `updated_paths_test_fetch.tsv`: Enriched version with local or `gs://` paths
- Local copies of all files in the specified output directory
- Optional: Uploaded files in GCS bucket

---
