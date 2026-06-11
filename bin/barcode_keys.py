def qualify_barcodes(barcodes, batch_key):
    """Return cell barcodes qualified with a stable cross-modality batch key."""
    suffix = str(batch_key)
    return [f"{barcode}_{suffix}" for barcode in map(str, barcodes)]
