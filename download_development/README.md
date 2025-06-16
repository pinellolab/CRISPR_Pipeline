# Generate_per_sample.py
```bash
python3  generate_per_sample.py --keypair igvf_key.json --accession IGVFDS9445RJOU --output test_fetch.tsv
```

```
The content of the .json file is :

{"key":"********","secret":"*********"}

```

#### Case the dataset does not have SeqSpec in the IGVF portal

```bash
python3 generate_per_sample.py --keypair igvf_key.json --accession IGVFDS7340YDHF --output test_fetch.tsv --hash_seqspec hash_seq_spec.yaml --rna_seqspec rna_seq_spec.yaml --sgrna_seqspec sgrna_seq_spec.yaml
```

                                        
# download_igvf.py
Use the previous extracted tsv file to download files and generate full paths


```bash
python3  download_igvf.py --sample test_fetch.tsv --access-key ***  --secret-key ***
```
 TODO:
 - [ ] download SeqSpecs    
