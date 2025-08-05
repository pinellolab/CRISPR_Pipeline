import csv
import requests
from requests.auth import HTTPBasicAuth
import os
import json
import argparse


PORTAL = 'https://api.data.igvf.org'

"""
Generate per-sample metadata TSV for Perturb-seq pipeline.

REQUIREMENTS:
-------------
1. Analysis Set:
    - Must specify all `input_file_sets` (measurement sets and auxiliary sets) to be processed together.
    - Must have exactly one `construct_library_sets`. This property is calculated from the underlying measurement sets and auxiliary sets.

2. Measurement Set:
    - Must include:
        - `strand_specificity`
        - `onlist_files`
        - `onlist_method` (must be "no combination")
    - Other `onlist_method` values are currently unsupported.

3. Auxiliary Sets:
    - Must include a `measurement_sets` link to its associated measurement set.
    - If `file_set_type` is `cell hashing barcode sequencing`, the `barcode_map`, i.e. the mapping between hashtags and barcodes is required.

4. Construct Library Set:
    - Must contain **exactly one** file in `integrated_content_files` with:
        - `content_type`: "guide RNA sequences"
        - `status`: one of ["in progress", "released", "preview"]
    - Deprecated guide files must be marked as "revoked", "archived", "deleted", or "replaced".

5. Sequence Files:
    - Must be of `content_type`: "reads"
    - Must have `illumina_read_type`: "R1" and "R2"
    - Must **not** have `status` in ["deleted", "revoked"]
    - `seqspecs`:
        - Must have a seqspec for each modality i.e. scRNA-seq, gRNA, hashing
        - seqspec must be linked to its sequence files through `seqspec_of`
        - seqspec `status` must be one of ["in progress", "released", "preview"]
        - seqspec must be validated i.e. `upload_status` = validated
    - All sequence files to be processed together should have the same seqspec read index per modality.
    - If `seqspecs` are missing, you must provide fallback YAML paths using CLI flags:
        --rna_seqspec, --sgrna_seqspec, --hash_seqspec

6. CLI Parameters:
    Required:
        --accession: Analysis Set accession (e.g. IGVFDS7340YDHF)
        --output:    Output path for generated `.tsv`

    Optional:
        --keypair:        Path to keypair JSON with `key` and `secret`
        --hash_seqspec:   Fallback YAML path for hashing
        --rna_seqspec:    Fallback YAML path for scRNA-seq
        --sgrna_seqspec:  Fallback YAML path for sgRNA
"""


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

    # Get guide RNA sequences reference
    guide_rna_sequences = ''
    construct_library_sets = analysis_set_object.get('construct_library_sets', [])
    if len(construct_library_sets) > 1:
        raise ValueError(
            f'Datasets with multiple guide libraries are not currently supported by the pipeline.'
        )
    construct_library_set_object = requests.get(f"{PORTAL}{construct_library_sets[0]}/@@embedded?format=json", auth=auth).json()
    integrated_content_files = construct_library_set_object.get('integrated_content_files', [])
    guide_rna_sequences = [
        file['accession']
        for file in integrated_content_files
        if (
            file['content_type'] == 'guide RNA sequences'
            and file.get('status') in ['in progress', 'preview', 'released']
        )
    ]    
    if len(guide_rna_sequences) > 1:
        raise ValueError(
            f'Datasets with multiple guide libraries are not currently supported by the pipeline.'
        )
    elif len(guide_rna_sequences) == 1:
        guide_rna_sequences = guide_rna_sequences[0]
    else:
        raise ValueError(
            f'Guide libraries are required for running this pipeline.'
        )
    
    measurement_sets_to_propertes = {}
    input_file_sets = sorted([file_set for file_set in analysis_set_object.get('input_file_sets', [])
                            if file_set.startswith('/measurement-sets/') or file_set.startswith('/auxiliary-sets/')],
                            reverse=True) # start with the measurement sets to create a mapping for measurement set-only properties
    for input_file_set in input_file_sets:
        print(f'Parsing input file set: {input_file_set}')
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
                raise ValueError(f'Datasets with onlist_method {onlist_method} are not currently supported by the pipeline.')
            strand_specificity = file_set_object.get('strand_specificity', '')
            
            if file_set_object['@id'] not in measurement_sets_to_propertes:
                measurement_sets_to_propertes[file_set_object['@id']] = {}
            measurement_sets_to_propertes[file_set_object['@id']]['barcode_onlist'] = barcode_onlist
            measurement_sets_to_propertes[file_set_object['@id']]['onlist_method'] = onlist_method
            measurement_sets_to_propertes[file_set_object['@id']]['strand_specificity'] = strand_specificity
        else:
            measurement_sets = file_set_object.get('measurement_sets', [])
            barcode_onlist = measurement_sets_to_propertes[measurement_sets[0]]['barcode_onlist']
            onlist_method = measurement_sets_to_propertes[measurement_sets[0]]['onlist_method']
            strand_specificity = measurement_sets_to_propertes[measurement_sets[0]]['strand_specificity']
            measurement_sets = ', '.join([measurement_set.split('/')[-2] for measurement_set in measurement_sets])


        # Get hashtag reference
        barcode_to_hashtag_map = ''
        if modality == 'cell hashing barcode sequencing':
            barcode_to_hashtag_map = file_set_object.get('barcode_map', '')
            if not(barcode_to_hashtag_map):
                raise ValueError(f'Missing barcode_map on auxiliary set {file_set_object["@id"]}. Datasets with cell hashing are required to specify a barcode to hash map.')
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
            seqspec_path = ''
            for seqspec in seqspecs:
                seqspec_file_object = requests.get(f"{PORTAL}{seqspec}/@@object?format=json", auth=auth).json()
                if seqspec_file_object.get('upload_status') != 'validated' or seqspec_file_object.get('status') not in ['in progress', 'preview', 'released']:
                    continue
                else:
                    seqspec_path = seqspec.split('/')[-2]
            if not seqspec_path:
                if hash_seqspec and rna_seqspec and sgrna_seqspec:
                    seqspec_path = modality_to_fallback_seqspec(
                        modality, hash_seqspec, rna_seqspec, sgrna_seqspec
                    )
                else:
                    raise ValueError(
                        f'Missing seqspec for modality {modality} (R1 {read1["@id"]}). '
                        'Provide a fallback using --hash_seqspec / --rna_seqspec / --sgrna_seqspec or make sure to deposit your seqspec in the data portal.'
                    )

            if not(strand_specificity):
                raise ValueError(f'Missing an associated strand specificty on {file_set_object["@id"]}. Please submit to the strand_specificity property on the data portal under the measurement set.')
            
            if not(barcode_onlist):
                raise ValueError(f'Missing an associated barcode onlist on {file_set_object["@id"]}. Please submit to the onlist_files property on the data portal under the measurement set.')

            if not(onlist_method):
                raise ValueError(f'Missing an associated onlist method on {file_set_object["@id"]}. Please submit to the onlist_method property on the data portal under the measurement set.')

            row = {
                'R1_path': read1['@id'].split('/')[-2],
                'R1_md5sum': read1['content_md5sum'],
                'R2_path': read2['@id'].split('/')[-2],
                'R2_md5sum': read2['content_md5sum'],
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
        raise ValueError("At least one row is missing a seqspec path – aborting.")

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
