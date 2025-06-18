import os
import requests
import pandas as pd
import argparse
import hashlib
import json
from tqdm import tqdm
from requests.auth import HTTPBasicAuth

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

def download_and_verify_file(accession: str, expected_md5: str, fastq_dir: str, auth: HTTPBasicAuth):
    """
    Downloads a single FASTQ file, shows progress, and verifies its MD5 checksum.
    Skips the download if the file already exists and the checksum is correct.

    Args:
        accession: The accession ID of the file to download.
        expected_md5: The expected MD5 checksum for verification.
        fastq_dir: The directory to save the file in.
        auth: HTTPBasicAuth object for authentication.
    """
    base_url = "https://api.data.igvf.org/sequence-files"
    download_url = f"{base_url}/{accession}/@@download/{accession}.fastq.gz"
    output_path = os.path.join(fastq_dir, f"{accession}.fastq.gz")

    # --- 1. Check if file already exists ---
    if os.path.exists(output_path):
        colored_print(Color.YELLOW, f"File {accession}.fastq.gz already exists. Verifying checksum...")
        # --- 3. Verify existing file ---
        actual_md5 = calculate_md5(output_path)
        if actual_md5 == expected_md5:
            colored_print(Color.GREEN, f"Checksum for existing file {accession} is correct. Skipping download.")
            return
        else:
            colored_print(Color.RED, f"Checksum mismatch for existing file {accession}. Expected {expected_md5}, got {actual_md5}. Re-downloading.")

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

    # --- 3. Verify the downloaded file ---
    colored_print(Color.YELLOW, f"Verifying checksum for {accession}...")
    actual_md5 = calculate_md5(output_path)
    if actual_md5 == expected_md5:
        colored_print(Color.GREEN, f"Checksum for {accession} is correct: {actual_md5}")
    else:
        colored_print(Color.RED, f"CRITICAL: Checksum mismatch for {accession}!")
        colored_print(Color.RED, f"  - Expected: {expected_md5}")
        colored_print(Color.RED, f"  - Got:      {actual_md5}")
        colored_print(Color.RED, f"The downloaded file at {output_path} may be corrupt.")


def process_sample_sheet(df, auth: HTTPBasicAuth) -> pd.DataFrame:
    """
    Processes a sample sheet to download all specified FASTQ file pairs and
    returns a new DataFrame with local file paths.

    Args:
        df: Pandas DataFrame with R1/R2 paths and md5sums.
        auth: HTTPBasicAuth object for authentication.

    Returns:
        A new Pandas DataFrame with R1_path and R2_path updated to full local paths.
    """
    fastq_dir = "fastq_files"
    os.makedirs(fastq_dir, exist_ok=True)

    local_paths_df = df.copy()

    r1_local_paths = []
    r2_local_paths = []

    # Use tqdm to create an overall progress bar for the file pairs
    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing all file pairs"):
        r1_accession = row['R1_path']
        r1_md5 = row['R1_md5sum']
        r2_accession = row['R2_path']
        r2_md5 = row['R2_md5sum']

        print("\n" + "-" * 60)
        colored_print(Color.GREEN, f"Processing pair: {r1_accession} & {r2_accession}")

        # Download and verify R1
        download_and_verify_file(r1_accession, r1_md5, fastq_dir, auth)

        # Download and verify R2
        download_and_verify_file(r2_accession, r2_md5, fastq_dir, auth)

        # Construct absolute paths for the new sample sheet
        r1_local_path = os.path.abspath(os.path.join(fastq_dir, f"{r1_accession}.fastq.gz"))
        r2_local_path = os.path.abspath(os.path.join(fastq_dir, f"{r2_accession}.fastq.gz"))
        r1_local_paths.append(r1_local_path)
        r2_local_paths.append(r2_local_path)

    # Update the DataFrame with the new local paths
    local_paths_df['R1_path'] = r1_local_paths
    local_paths_df['R2_path'] = r2_local_paths

    return local_paths_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Download and verify FASTQ files from the IGVF portal.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example usage:
# Using keypair JSON file:
python3 download_igvf.py --sample test_fetch.tsv --keypair keypair.json
"""
    )
    parser.add_argument('--sample', required=True,
                        help='Path to the TSV file containing sample information.')
    parser.add_argument('--keypair',
                        help='Path to JSON file containing IGVF-API key pair.')

    args = parser.parse_args()

    try:
        # Get authentication using the get_auth function
        auth = get_auth(args.keypair)

        df = pd.read_csv(args.sample, sep='\t')
        required_cols = ['R1_path', 'R2_path', 'R1_md5sum', 'R2_md5sum']

        # Check for required columns
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"TSV file is missing required columns: {', '.join(missing_cols)}")

        # Start the download process and get the DataFrame with local paths
        local_paths_df = process_sample_sheet(df, auth)

        # --- Create the new TSV file with local paths ---
        output_basename = os.path.basename(args.sample)
        local_paths_filename = f"local_paths_{output_basename}"

        local_paths_df.to_csv(local_paths_filename, sep='\t', index=False)

        colored_print(Color.GREEN, "\nAll download and verification tasks complete.")
        colored_print(Color.GREEN, f"Successfully created sample sheet with local paths: {local_paths_filename}")

    except FileNotFoundError:
        colored_print(Color.RED, f"Error: The file '{args.sample}' was not found.")
    except pd.errors.EmptyDataError:
        colored_print(Color.RED, f"Error: The file '{args.sample}' is empty or improperly formatted.")
    except RuntimeError as e:
        colored_print(Color.RED, f"Authentication error: {str(e)}")
    except Exception as e:
        colored_print(Color.RED, f"An unexpected error occurred: {str(e)}")
