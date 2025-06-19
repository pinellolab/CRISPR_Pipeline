
process seqSpecParser {
    cache 'lenient'

    input:
    path seqSpec_yaml
    path barcode_file
    val modalities

    output:
    path "${modalities}_parsed_seqSpec.txt", emit: parsed_seqspec
    path "barcode_file.txt", emit: barcode_file

    script:
    """
    parsing_guide_metadata.py --modalities ${modalities} --yaml_file ${seqSpec_yaml} --whitelist ${barcode_file} --output_file ${modalities}_parsed_seqSpec.txt
    echo "Path to the barcode file: ${barcode_file}"
    cp "${barcode_file}" barcode_file.txt
    """
}
