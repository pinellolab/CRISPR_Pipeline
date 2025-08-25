import os
import requests
import pandas as pd
import argparse
import hashlib
import json
import gzip
import shutil
from tqdm import tqdm
from requests.auth import HTTPBasicAuth
from google.cloud import storage


SEQUENCE_FILE_COLUMNS = ['R1_path', 'R2_path']
CONFIGURATION_FILE_COLUMNS = ['seqspec']
TABULAR_FILE_COLUMNS = ['barcode_onlist', 'guide_design', 'barcode_hashtag_map']

# Helper class for colored console output
class Color:
    """A simple class to add color to terminal output."""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'

def get_auth(keypair_path=None):
    if keypair_path:
        with open(keypair_path) as f:
            keypair = json.load(f)
            return HTTPBasicAuth(keypair['key'], keypair['secret'])

    key = os.getenv('IGVF_API_KEY')
    secret = os.getenv('IGVF_SECRET_KEY')
    if key and secret:
        return HTTPBasicAuth(key, secret)

    raise RuntimeError('No credentials provided. Set IGVF_API_KEY and IGVF_SECRET_KEY or provide a keypair JSON.')

def colored_print(color, message):
    """Prints a message in the specified color."""
    print(f"{color}{message}{Color.END}")

def calculate_md5(filepath: str, chunk_size: int = 8192) -> str:
    """
    Calculate the MD5 checksum of a file.
    """
    md5 = hashlib.md5()
    try:
        with open(filepath, 'rb') as f:
            while chunk := f.read(chunk_size):
                md5.update(chunk)
    except IOError as e:
        colored_print(Color.RED, f"Error reading file for MD5 calculation: {e}")
        return ""
    return md5.hexdigest()

def download_and_verify_file(accession: str, expected_md5: str, download_dir: str,
                             auth: HTTPBasicAuth, file_type: str, args):
    """
    Downloads one file from IGVF, optionally gunzips, verifies checksum (if provided),
    and optionally uploads the resulting artifact to GCS.

    Returns:
        (local_target_path, gcs_url_or_None)
    """
    # Determine download URL and file extension based on file type
    file_type_to_url = {
        'sequence': ('sequence-files', 'fastq.gz'),
        'tabular': ('tabular-files', 'tsv.gz'),
        'configuration': ('configuration-files', 'yaml.gz')
    }
    if file_type not in file_type_to_url:
        raise ValueError(f"Unknown file_type '{file_type}'")


    route, ext = file_type_to_url[file_type]
    base_url = f"https://api.data.igvf.org/{route}"
    download_url = f"{base_url}/{accession}/@@download/{accession}.{ext}"

    # Paths
    output_gz_path = os.path.join(download_dir, f"{accession}.{ext}")
    unzipped_name = f"{accession}.{ext[:-3]}" if ext.endswith('.gz') else f"{accession}.{ext}"
    unzipped_path = os.path.join(download_dir, unzipped_name)

    # Decide which artifact we want to end up with
    want_unzipped = bool(args.gunzip)
    target_path = unzipped_path if want_unzipped else output_gz_path
    gcs_url = None

    if args.dry_run:
        colored_print(Color.YELLOW, f"[Dry run] Would download: {download_url}")
        colored_print(Color.YELLOW, f"[Dry run] Would produce: {target_path} "
                                    f"({'unzipped' if want_unzipped else 'gzipped'})")
        return None, None

    # If target already exists, verify (if requested) and (optionally) upload
    if os.path.exists(target_path):
        if expected_md5:
            colored_print(Color.YELLOW, f"Verifying checksum for existing {accession} "
                                        f"({'unzipped' if want_unzipped else 'gzipped'})...")
            actual_md5 = calculate_md5(target_path)
            if actual_md5 != expected_md5:
                colored_print(Color.RED, "Checksum mismatch on existing file... Re-downloading...")
            else:
                colored_print(Color.GREEN, f"Checksum is correct. Skipping download for {accession}.")
                if args.gcs_upload:
                    gcs_url = upload_to_gcs(target_path, args)
                return target_path, gcs_url
        else:
            colored_print(Color.YELLOW, f"No checksum provided for {accession}. Using existing file.")
            if args.gcs_upload:
                gcs_url = upload_to_gcs(target_path, args)
            return target_path, gcs_url

    # --- Download gz payload first (API serves .gz) ---
    try:
        response = requests.get(download_url, auth=auth, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with open(output_gz_path, 'wb') as f, tqdm(
            desc=f"Downloading {accession}",
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
            leave=False
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    size = f.write(chunk)
                    bar.update(size)

        if total_size != 0 and os.path.getsize(output_gz_path) != total_size:
            raise IOError("Download incomplete: Size mismatch.")
    except requests.exceptions.RequestException as e:
        colored_print(Color.RED, f"Failed to download {accession}: {e}")
        if os.path.exists(output_gz_path):
            os.remove(output_gz_path)
        return None, None
    except IOError as e:
        colored_print(Color.RED, f"Error writing file {accession}: {e}")
        if os.path.exists(output_gz_path):
            os.remove(output_gz_path)
        return None, None

    colored_print(Color.GREEN, f"Downloaded {accession} to {os.path.basename(output_gz_path)}")

    # --- Optionally gunzip ---
    produced_path = output_gz_path
    if want_unzipped:
        try:
            with gzip.open(output_gz_path, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(output_gz_path)
            produced_path = unzipped_path
            colored_print(Color.GREEN, f"Gunzipped {accession} to {os.path.basename(unzipped_path)}")
        except Exception as e:
            colored_print(Color.RED, f"Failed to gunzip {accession}: {e}")
            produced_path = output_gz_path  # fall back to gz

    # --- Verify checksum against produced artifact (if provided) ---
    if expected_md5:
        artifact_label = 'unzipped' if produced_path == unzipped_path else 'gzipped'
        colored_print(Color.YELLOW, f"Verifying checksum for {accession} ({artifact_label})...")
        actual_md5 = calculate_md5(produced_path)
        if actual_md5 == expected_md5:
            colored_print(Color.GREEN, f"Checksum for {accession} is correct: {actual_md5}")
        else:
            colored_print(Color.RED, f"CRITICAL: Checksum mismatch for {accession}!")
            colored_print(Color.RED, f"  - Expected: {expected_md5}")
            colored_print(Color.RED, f"  - Got:      {actual_md5}")
            colored_print(Color.RED, f"The downloaded file at {produced_path} may be corrupt.")
    else:
        colored_print(Color.YELLOW, f"No checksum provided for {accession}. Skipping verification.")

    # --- Upload to GCS if enabled ---
    gcs_url = upload_to_gcs(produced_path, args) if args.gcs_upload else None
    return produced_path, gcs_url

def process_sample_sheet(df, auth: HTTPBasicAuth, args, file_types='all', output_dir='downloads') -> pd.DataFrame:
    """
    For each IGVFFI accession in columns, download (if needed), verify (if MD5 provided),
    gunzip (if --gunzip), and optionally upload to GCS. Then update the original columns in-place with the final path:
      - if --gcs-upload: gs://bucket/prefix/filename
      - else: local filesystem path

    Columns handled:
      - FASTQ: R1_path, R2_path (when file_types in ['fastq','all'])
      - Configuration: seqspec (when file_types in ['other','all'])
      - Tabular: barcode_onlist, guide_design, barcode_hashtag_map (when file_types in ['other','all'])
    """
    os.makedirs(output_dir, exist_ok=True)

    cols_to_update = []
    if file_types in ['fastq', 'all']:
        cols_to_update.extend(SEQUENCE_FILE_COLUMNS)
    if file_types in ['other', 'all']:
        cols_to_update.extend(CONFIGURATION_FILE_COLUMNS)
        cols_to_update.extend(TABULAR_FILE_COLUMNS)

    col_to_file_type = {}
    for col in SEQUENCE_FILE_COLUMNS:
        col_to_file_type[col] = 'sequence'
    for col in CONFIGURATION_FILE_COLUMNS:
        col_to_file_type[col] = 'configuration'
    for col in TABULAR_FILE_COLUMNS:
        col_to_file_type[col] = 'tabular'

    # Cache: accession -> resolved path (local or gs://) so we don't re-download
    resolved_path_cache = {}

    # Walk rows and update in-place
    for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc="Downloading/updating paths"):
        for col in cols_to_update:
            if col not in df.columns:
                continue
            value = row[col]

            # Only touch IGVF file accessions; leave non-IGVFFI values (already paths) as-is
            if not (pd.notna(value) and isinstance(value, str) and value.startswith("IGVFFI")):
                continue

            accession = value
            # If already processed, reuse the resolved path
            if accession in resolved_path_cache:
                df.at[idx, col] = resolved_path_cache[accession]
                continue

            # For FASTQs, try to pick up MD5 if present
            expected_md5 = None
            if col in SEQUENCE_FILE_COLUMNS:
                md5_col = f"{col.replace('_path', '')}_md5sum"
                if md5_col in df.columns:
                    expected_md5 = row.get(md5_col)

            file_type = col_to_file_type[col]
            result = download_and_verify_file(
                accession,
                expected_md5,
                output_dir,
                auth,
                file_type=file_type,
                args=args
            )

            # In dry-run or error cases, keep original value
            if not result:
                continue

            local_path, gcs_path = result
            final_path = gcs_path if (args.gcs_upload and gcs_path) else local_path

            # If failed to resolve a final path, don't clobber the accession
            if final_path:
                resolved_path_cache[accession] = final_path
                df.at[idx, col] = final_path

    return df

def upload_to_gcs(path, args):
    """Upload file to GCS if it doesn't already exist or if forced.
    Returns the full gs:// URL of the uploaded (or existing) object, or None on failure.
    """
    gcs_rel_path = os.path.join(args.gcs_prefix, os.path.basename(path))
    client = storage.Client()
    bucket = client.bucket(args.gcs_bucket)
    blob = bucket.blob(gcs_rel_path)

    full_uri = f"gs://{args.gcs_bucket}/{gcs_rel_path}"

    if blob.exists() and not args.gcs_force_upload:
        print(f"Skipping GCS upload (already exists): {full_uri}")
        return full_uri

    print(f"Uploading to GCS: {full_uri}")
    try:
        blob.upload_from_filename(path)
        return full_uri
    except Exception as e:
        print(f"Failed to upload {os.path.basename(path)} to GCS: {e}")
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Download and verify FASTQ, configuration, and tabular files from the IGVF portal. '
                    'Also, optionally upload the files to a google cloud bucket.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
        Example usage:

        # Download everything (FASTQs, seqspecs, etc.), gunzip, and upload to GCS
        python3 download_igvf.py --sample per_sample.tsv --keypair keypair.json --gunzip --gcs-upload --gcs-bucket igvf-bkt --gcs-prefix pipeline_inputs

        # Download only FASTQ files and gunzip
        python3 download_igvf.py --sample per_sample.tsv --keypair keypair.json --file-types fastq --gunzip --output-dir ./fastq_dir

        # Download only non-FASTQ files, keep compressed (default)
        python3 download_igvf.py --sample per_sample.tsv --keypair keypair.json --file-types other --output-dir ./reference_files
        """
    )
    parser.add_argument('--sample', required=True,
                        help='Path to the TSV file containing sample information.')
    parser.add_argument('--keypair',
                        help='Path to JSON file containing IGVF-API key pair.')
    parser.add_argument('--file-types', choices=['fastq', 'other', 'all'], default='all',
                        help='Which file types to download: fastq, others (including tsvs and yaml files), or all (default).')
    parser.add_argument('--output-dir', default='downloads',
                        help='Directory to store downloaded files (default: ./downloads)')

    parser.add_argument('--gunzip', action='store_true', default=False,
                        help='Gunzip downloaded files (default: keep compressed).')

    parser.add_argument('--gcs-upload', action='store_true',
                        help='Upload downloaded files to GCS. Requires --gcs-bucket and --gcs-prefix. '
                             'If specified, the updated paths in the sample sheet will include gs:// URLs.')
    parser.add_argument('--gcs-bucket',
                        help='Name of the GCS bucket to upload to.')
    parser.add_argument('--gcs-prefix',
                        help='(Required if --gcs-upload) Subfolder inside the GCS bucket '
                             '(e.g., pipeline_runs/IGVFDS6673ZFFG_2025_08_06).')
    parser.add_argument('--gcs-force-upload', action='store_true',
                        help='Force upload to GCS even if file already exists in the bucket')

    parser.add_argument('--dry-run', action='store_true',
                        help='Print actions without performing download/upload.')

    args = parser.parse_args()
    if args.gcs_upload and (not args.gcs_bucket or not args.gcs_prefix):
        raise ValueError("When using --gcs-upload, both --gcs-bucket and --gcs-prefix must be provided.")

    try:
        # Get authentication using the get_auth function
        auth = get_auth(args.keypair)

        df = pd.read_csv(args.sample, sep='\t')

        # Only require FASTQ-related columns if downloading FASTQs
        if args.file_types in ['fastq', 'all']:
            required_cols = ['R1_path', 'R2_path', 'R1_md5sum', 'R2_md5sum']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                raise ValueError(f"TSV file is missing required FASTQ columns: {', '.join(missing_cols)}")

        # Start the download process and get the DataFrame with local/GCS paths
        local_paths_df = process_sample_sheet(df, auth, args, file_types=args.file_types, output_dir=args.output_dir)

        # Create a new TSV with additional columns for local/GCS file paths
        output_basename = os.path.basename(args.sample)
        updated_paths_filename = f"updated_paths_{output_basename}"
        local_paths_df.to_csv(updated_paths_filename, sep='\t', index=False)

        colored_print(Color.GREEN, "\nAll download and verification tasks complete.")
        colored_print(Color.GREEN, f"Successfully created sample sheet with updated paths: {updated_paths_filename}")

        if args.gcs_upload:
            try:
                gcs_url = upload_to_gcs(updated_paths_filename, args)
                if gcs_url:
                    colored_print(Color.GREEN, f"Sample sheet uploaded to: {gcs_url}")
                else:
                    colored_print(Color.RED, "Sample sheet upload to GCS failed.")
            except Exception as e:
                colored_print(Color.RED, f"Failed to upload sample sheet to GCS: {e}")

    except FileNotFoundError as e:
        colored_print(Color.RED, f"[FileNotFoundError] {e}")
        raise
    except pd.errors.EmptyDataError as e:
        colored_print(Color.RED, f"[EmptyDataError] {e}")
        raise
    except RuntimeError as e:
        colored_print(Color.RED, f"Authentication error: {str(e)}")
        raise
    except Exception as e:
        colored_print(Color.RED, f"[Unexpected {type(e).__name__}] {e}")
        raise
