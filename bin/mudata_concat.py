#!/usr/bin/env python
import argparse
import os
import mudata as md



def filter_genes_by_cells(mdata, min_cells_faction):
    """
        Filter genes in the MuData object based on the minimum number of cells they must be expressed in.
    """
    
    index_filter = ((mdata['gene'].X >0).sum(0).A1 > int(mdata['gene'].n_obs * min_cells_faction)  ) #gene with more than 10k cells
    mdata.mod['gene'] = mdata.mod['gene'][:, index_filter ]
    return  mdata

def concat_mudatas(input_files, output_file, min_cells_fraction=0.05):
    """
    Concatenate multiple MuData files. If only one file is provided, it's copied to the output.
    """
    files = sorted(input_files, key=lambda x: os.path.basename(x))
    if not files:
        print(f"No files found: {input_files}")
        return

    print(f"Found {len(files)} files to concatenate")

    # Handle single file case
    if len(files) == 1:
        print(f"Only one file found. Copying {files[0]} to {output_file}")
        single_mdata = md.read(files[0])
        single_mdata = filter_genes_by_cells(single_mdata, min_cells_fraction)  # Filter genes based on minimum cells
        print(f"Saving MuData with {single_mdata.n_obs} cells to {output_file}")
        single_mdata.write(output_file)
        return

    # Handle multiple files case
    print("Concatenating all MuData objects...")
    combined_mdata = md.concat([md.read(ff) for ff in files], merge='first', uns_merge='first', join='outer')


    print ('filtering genes')
    combined_mdata = filter_genes_by_cells(combined_mdata, min_cells_fraction)  # Filter genes based on minimum cells


    print(f"Saving combined MuData with {combined_mdata.n_obs} cells to {output_file}")
    combined_mdata.write(output_file)

    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Concatenate MuData files")
    parser.add_argument("-i", "--input", dest="input", nargs="+", required=True, help="Input mudata files")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output file path")
    parser.add_argument("-g", "--gene_filter", dest="filter genes based in a percentage", type=float, default=0.05, help="Minimum number of cells a gene must be expressed in to be retained")
    args = parser.parse_args()

    concat_mudatas(args.input, args.output, args.gene_filter)

if __name__ == "__main__":
    main()
