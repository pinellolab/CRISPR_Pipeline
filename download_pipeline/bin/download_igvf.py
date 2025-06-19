import os
import requests
import pandas as pd
import argparse
import hashlib
import json
from tqdm import tqdm
from requests.auth import HTTPBasicAuth


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

    Args:
        filepath: Path to the file.
        chunk_size: The size of chunks to read from the file.

    Returns:
        The MD5 checksum in hexadecimal format.
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

def download_and_verify_file(accession: str, expected_md5: str, download_dir: str, auth: HTTPBasicAuth, file_type: str):
    """
        Downloads a single file (FASTQ, tabular, or configuration), shows progress,
        and verifies its MD5 checksum if provided. Skips the download if the file
        already exists and is verified (or if no checksum is provided).

    Args:
        accession: The accession ID of the file to download.
        expected_md5: The expected MD5 checksum for verification.
        dir: The directory to save the file in.
        auth: HTTPBasicAuth object for authentication.
        file_type: file object type being downloaded (i.e. sequence, tabular, or configuration)
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
    output_path = os.path.join(download_dir, f"{accession}.{ext}")

    # --- 1. Check if file already exists ---
    if os.path.exists(output_path):
        if expected_md5 is None:
            colored_print(Color.YELLOW, f"No checksum provided for {accession}. Skipping verification and download.")
            return
        else:
            colored_print(Color.YELLOW, f"Verifying checksum for {accession}...")
            actual_md5 = calculate_md5(output_path)
            if actual_md5 == expected_md5:
                colored_print(Color.GREEN, f"Checksum for existing file {accession} is correct. Skipping download.")
                return
            else:
                colored_print(Color.RED, f"CRITICAL: Checksum mismatch for existing file {accession}!")
                colored_print(Color.RED, f"  - Expected: {expected_md5}")
                colored_print(Color.RED, f"  - Got:      {actual_md5}")
                colored_print(Color.RED, f"Re-downloading {accession}...")

    # --- 2. Download the file with a progress bar ---
    try:
        response = requests.get(download_url, auth=auth, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        total_size = int(response.headers.get('content-length', 0))

        with open(output_path, 'wb') as f, tqdm(
            desc=f"Downloading {accession}",
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
            leave=False # Prevents this progress bar from staying after completion
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                size = f.write(chunk)
                bar.update(size)

        if total_size != 0 and bar.n != total_size:
            raise IOError("Download incomplete: Size mismatch.")

    except requests.exceptions.RequestException as e:
        colored_print(Color.RED, f"Failed to download {accession}: {e}")
        if os.path.exists(output_path):
            os.remove(output_path) # Clean up partial download
        return
    except IOError as e:
        colored_print(Color.RED, f"Error writing file {accession}: {e}")
        if os.path.exists(output_path):
            os.remove(output_path) # Clean up partial download
        return

    colored_print(Color.GREEN, f"Downloaded {accession} successfully.")

    # --- 3. Verify the downloaded file if expected_md5 is provided ---
    if expected_md5 is not None:
        colored_print(Color.YELLOW, f"Verifying checksum for {accession}...")
        actual_md5 = calculate_md5(output_path)
        if actual_md5 == expected_md5:
            colored_print(Color.GREEN, f"Checksum for {accession} is correct: {actual_md5}")
        else:
            colored_print(Color.RED, f"CRITICAL: Checksum mismatch for {accession}!")
            colored_print(Color.RED, f"  - Expected: {expected_md5}")
            colored_print(Color.RED, f"  - Got:      {actual_md5}")
            colored_print(Color.RED, f"The downloaded file at {output_path} may be corrupt.")
    else:
        colored_print(Color.YELLOW, f"No checksum provided for {accession}. Skipping verification.")


def process_sample_sheet(df, auth: HTTPBasicAuth, file_types='all', output_dir='downloads') -> pd.DataFrame:

    os.makedirs(output_dir, exist_ok=True)
    local_paths_df = df.copy()

    # Store local file paths
    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Downloading files per sample"):
        for col in df.columns:
            if file_types in ['fastq', 'all'] and col in SEQUENCE_FILE_COLUMNS:
                accession = row[col]
                md5_col = f"{col}_md5sum"
                expected_md5 = row.get(md5_col)
                if pd.notna(accession) and isinstance(accession, str) and accession.startswith("IGVFFI"):
                    download_and_verify_file(accession, expected_md5, output_dir, auth, file_type='sequence')
                    local_paths_df.at[_, f'local_{col}'] = os.path.abspath(os.path.join(output_dir, f"{accession}.fastq.gz"))

            elif file_types in ['other', 'all']:
                if col in CONFIGURATION_FILE_COLUMNS and pd.notna(row[col]):
                    accession = row[col]
                    download_and_verify_file(accession, None, output_dir, auth, file_type='configuration')
                    local_paths_df.at[_, f'local_{col}'] = os.path.abspath(os.path.join(output_dir, f"{accession}.yaml.gz"))

                elif col in TABULAR_FILE_COLUMNS and pd.notna(row[col]):
                    accession = row[col]
                    download_and_verify_file(accession, None, output_dir, auth, file_type='tabular')
                    local_paths_df.at[_, f'local_{col}'] = os.path.abspath(os.path.join(output_dir, f"{accession}.tsv.gz"))

    return local_paths_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Download and verify FASTQ, configuration, and tabular files from the IGVF portal.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
        Example usage:

        # Download everything (FASTQs, seqspecs, barcode maps, guide designs, etc.)
        python3 download_igvf.py --sample per_sample.tsv --keypair keypair.json --output-dir ./pipeline_inputs

        # Download only FASTQ files (R1/R2)
        python3 download_igvf.py --sample per_sample.tsv --keypair keypair.json --file-types fastq --output-dir ./fastq_dir

        # Download only non-FASTQ files (seqspec, barcode_onlist, guide_design, barcode_hashtag_map)
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

    args = parser.parse_args()


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

        # Start the download process and get the DataFrame with local paths
        local_paths_df = process_sample_sheet(df, auth, file_types=args.file_types, output_dir=args.output_dir)

        # --- Create a new TSV with additional columns for local file paths ---
        output_basename = os.path.basename(args.sample)
        local_paths_filename = f"local_paths_{output_basename}"

        local_paths_df.to_csv(local_paths_filename, sep='\t', index=False)

        colored_print(Color.GREEN, "\nAll download and verification tasks complete.")
        colored_print(Color.GREEN, f"Successfully created sample sheet with local paths: {local_paths_filename}")

    except FileNotFoundError as e:
        colored_print(Color.RED, f"[FileNotFoundError] {e}")
        raise
    except pd.errors.EmptyDataError as e:
        colored_print(Color.RED, f"[EmptyDataError] {e}")
        raise
    except Exception as e:
        colored_print(Color.RED, f"[Unexpected {type(e).__name__}] {e}")
        raise
    except RuntimeError as e:
        colored_print(Color.RED, f"Authentication error: {str(e)}")
        raise
