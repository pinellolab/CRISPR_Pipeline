import csv
import requests
from requests.auth import HTTPBasicAuth
import os
import json
import argparse


PORTAL = 'https://api.data.igvf.org'


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


def write_tsv(data_list, output_path):
    if not data_list:
        raise ValueError("No data provided for TSV output.")
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=data_list[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(data_list)


def modality_to_fallback_seqspec(modality, hash_seqspec, rna_seqspec, sgrna_seqspec):
    modality = modality.lower().strip()
    if modality == 'scrna sequencing':
        return rna_seqspec
    elif modality == 'grna sequencing':
            return sgrna_seqspec
    elif modality == 'cell hashing barcode sequencing':
            return hash_seqspec
    return None


def generate_per_sample_tsv(analysis_set_accession, output_path, auth, hash_seqspec=None, rna_seqspec=None, sgrna_seqspec=None):
    response = requests.get(
        f'{PORTAL}/analysis-sets/{analysis_set_accession}/@@object?format=json', auth=auth
    )
    analysis_set_object = response.json()
    per_sample_rows = []

    for input_file_set in analysis_set_object.get('input_file_sets', []):
        if input_file_set.startswith('/measurement-sets/') or input_file_set.startswith('/auxiliary-sets/'):
            file_set_object = requests.get(f"{PORTAL}{input_file_set}/@@object?format=json", auth=auth).json()
            file_set_type = file_set_object['file_set_type']
            modality = 'scRNA sequencing' if file_set_type == 'experimental data' else file_set_type
            accession = file_set_object['accession']
            files = file_set_object.get('files', [])

            barcode_onlist = ''
            strand_specificity = ''
            if file_set_object['@id'].startswith('/measurement-sets/'):
                measurement_sets = accession
                barcode_onlist = file_set_object.get('onlist_files', [])
                onlist_method = file_set_object.get('onlist_method', '')
                if onlist_method and onlist_method != 'no combination':
                    print(f'{onlist_method} not currently supported') # Upgrade this to a ValueError later
                strand_specificity = file_set_object.get('strand_specificity', '')
            else:
                measurement_sets = ', '.join([measurement_set.split('/')[-2] for measurement_set in file_set_object.get('measurement_sets', [])])

            # Get guide RNA sequences reference
            guide_rna_sequences = ''
            if modality == 'gRNA sequencing':
                construct_library_sets = file_set_object.get('construct_library_sets', [])
                if len(construct_library_sets) > 1:
                    raise ValueError(
                        f'Datasets with multiple guide libraries are not currently supported by the pipeline.'
                    )
                construct_library_set_object = requests.get(f"{PORTAL}{construct_library_sets[0]}/@@embedded?format=json", auth=auth).json()
                integrated_content_files = construct_library_set_object.get('integrated_content_files', [])
                guide_rna_sequences = [integrated_content_file['accession'] for integrated_content_file in integrated_content_files if integrated_content_file['content_type'] == 'guide RNA sequences']
                if len(guide_rna_sequences) > 1:
                    raise ValueError(
                        f'Datasets with multiple guide libraries are not currently supported by the pipeline.'
                    )
                elif len(guide_rna_sequences) == 1:
                    guide_rna_sequences = guide_rna_sequences[0]
                else:
                    raise ValueError(
                        f'Datasets without guide libraries are not supported by the pipeline.'
                    )


            # Get hashtag reference
            barcode_to_hashtag_map = ''
            if modality == 'cell hashing barcode sequencing':
                barcode_to_hashtag_map = file_set_object.get('barcode_map', '')
                if not(barcode_to_hashtag_map):
                    print('No barcode to hashtag mapping provided.') # Upgrade this to a ValueError later
                else:
                    barcode_to_hashtag_map = barcode_to_hashtag_map.split('/')[-2]

            sequence_file_index = {}
            for file in files:
                file_object = requests.get(f"{PORTAL}{file}/@@object?format=json", auth=auth).json()
                if not file_object['@id'].startswith('/sequence-files/') or file_object.get('status') in ['deleted', 'revoked'] or file_object['content_type'] != 'reads' or file_object.get('illumina_read_type') not in ['R1', 'R2']:
                    continue

                key = (
                    file_object.get('sequencing_run'),
                    file_object.get('lane'),
                    file_object.get('flowcell_id'),
                    file_object.get('index')
                )
                read_type = file_object.get('illumina_read_type')

                if key not in sequence_file_index:
                    sequence_file_index[key] = {}

                sequence_file_index[key][read_type] = file_object

            for key, reads in sequence_file_index.items():
                read1 = reads.get('R1')
                read2 = reads.get('R2')
                if not read1 or not read2:
                    continue

                # Handle seqspec fallback
                seqspecs = read1.get('seqspecs', [])
                if seqspecs and seqspecs[0]:
                    seqspec_path = seqspecs[0].split('/')[-2]
                else:
                    seqspec_path = modality_to_fallback_seqspec(
                        modality, hash_seqspec, rna_seqspec, sgrna_seqspec
                    )

                if not seqspec_path:
                    raise RuntimeError(
                        f"Missing seqspec for modality '{modality}' (R1 {read1['@id']}). "
                        "Provide a fallback using --hash_seqspec / --rna_seqspec / --sgrna_seqspec or make sure to deposited your seqspec in the portal"
                    )
                row = {
                    'R1_path': read1['@id'].split('/')[-2],
                    'R1_md5sum': read1['md5sum'],
                    'R2_path': read2['@id'].split('/')[-2],
                    'R2_md5sum': read2['md5sum'],
                    'file_modality': modality,
                    'file_set': accession,
                    'measurement_sets': measurement_sets,
                    'sequencing_run': key[0],
                    'lane': key[1],
                    'flowcell_id': key[2],
                    'index': key[3],
                    'seqspec': seqspec_path,
                    'barcode_onlist': barcode_onlist[0].split('/')[-2] if barcode_onlist else '',
                    'onlist_method': onlist_method,
                    'strand_specificity': strand_specificity,
                    'guide_design': guide_rna_sequences,
                    'barcode_hashtag_map': barcode_to_hashtag_map
                }

                per_sample_rows.append(row)

    # Final validation
    if any(not row['seqspec'] for row in per_sample_rows):
        raise RuntimeError("At least one row is missing a seqspec path – aborting.")

    write_tsv(per_sample_rows, output_path)


def main():
    parser = argparse.ArgumentParser(description="Generate per-sample TSV for Perturb-seq pipeline. Example usage:\n"
                                                 "python3 generate_per_sample.py --keypair igvf_key.json --accession IGVFDS7340YDHF "
                                                 "--output test_fetch.tsv --hash_seqspec hash_seq_spec.yaml "
                                                 "--rna_seqspec rna_seq_spec.yaml --sgrna_seqspec sgrna_seq_spec.yaml")
    parser.add_argument("--accession", required=True, help="Analysis set accession")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    parser.add_argument("--keypair", help="Optional path to keypair JSON file")
    parser.add_argument("--hash_seqspec", help="Fallback path to hash seqspec")
    parser.add_argument("--rna_seqspec", help="Fallback path to RNA seqspec")
    parser.add_argument("--sgrna_seqspec", help="Fallback path to sgRNA seqspec")
    args = parser.parse_args()

    print(f"Generating per-sample TSV from analysis set '{args.accession}'...")
    auth = get_auth(args.keypair)
    generate_per_sample_tsv(
        args.accession,
        args.output,
        auth,
        hash_seqspec=args.hash_seqspec,
        rna_seqspec=args.rna_seqspec,
        sgrna_seqspec=args.sgrna_seqspec
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        raise


#seqspec in the analysis set
#python3  generate_per_sample.py --keypair igvf_key.json --accession IGVFDS9445RJOU --output test_fetch.tsv
#seqspec not in the analysis set, but provided as fallback
#python3 generate_per_sample.py --keypair igvf_key.json --accession IGVFDS7340YDHF --output test_fetch.tsv --hash_seqspec hash_seq_spec.yaml --rna_seqspec rna_seq_spec.yaml --sgrna_seqspec sgrna_seq_spec.yaml
