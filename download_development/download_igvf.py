import os
import requests
import pandas as pd
import argparse
import hashlib
from tqdm import tqdm

# Helper class for colored console output
class Color:
    """A simple class to add color to terminal output."""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'

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

def download_and_verify_file(accession: str, expected_md5: str, fastq_dir: str, auth: tuple):
    """
    Downloads a single FASTQ file, shows progress, and verifies its MD5 checksum.
    Skips the download if the file already exists.

    Args:
        accession: The accession ID of the file to download.
        expected_md5: The expected MD5 checksum for verification.
        fastq_dir: The directory to save the file in.
        auth: A tuple containing the (access_key, secret_key) for authentication.
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


def process_sample_sheet(df, access_key: str, secret_key: str):
    """
    Processes a sample sheet to download all specified FASTQ file pairs.

    Args:
        df: Pandas DataFrame with R1/R2 paths and md5sums.
        access_key: IGVF API access key.
        secret_key: IGVF API secret key.
    """
    fastq_dir = "fastq_files"
    os.makedirs(fastq_dir, exist_ok=True)
    auth = (access_key, secret_key)

    # Use tqdm to create an overall progress bar for the file pairs
    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing all file pairs"):
        r1_accession = row['R1_path']
        r1_md5 = row['R1_md5sum']
        r2_accession = row['R2_path']
        r2_md5 = row['R2_md5sum']

        print("-" * 50)
        colored_print(Color.GREEN, f"Processing pair: {r1_accession} & {r2_accession}")
        
        # Download and verify R1
        download_and_verify_file(r1_accession, r1_md5, fastq_dir, auth)
        
        # Download and verify R2
        download_and_verify_file(r2_accession, r2_md5, fastq_dir, auth)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Download and verify FASTQ files from the IGVF portal.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example usage:
  python3 download_igvf.py --sample test_fetch.tsv --access-key YOUR_KEY --secret-key YOUR_SECRET
"""
    )
    parser.add_argument('--sample', required=True, 
                        help='Path to the TSV file containing sample information.')
    parser.add_argument('--access-key', required=True,
                        help='IGVF API access key.')
    parser.add_argument('--secret-key', required=True,
                        help='IGVF API secret key.')
    
    args = parser.parse_args()
    
    try:
        df = pd.read_csv(args.sample, sep='\t')
        required_cols = ['R1_path', 'R2_path', 'R1_md5sum', 'R2_md5sum']
        
        # Check for required columns
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"TSV file is missing required columns: {', '.join(missing_cols)}")
        
        # Start the download process
        process_sample_sheet(df, args.access_key, args.secret_key)
        
        colored_print(Color.GREEN, "\nAll download and verification tasks complete.")
        
    except FileNotFoundError:
        colored_print(Color.RED, f"Error: The file '{args.sample}' was not found.")
    except pd.errors.EmptyDataError:
        colored_print(Color.RED, f"Error: The file '{args.sample}' is empty or improperly formatted.")
    except Exception as e:
        colored_print(Color.RED, f"An unexpected error occurred: {str(e)}")