#!/bin/usr python
import pandas as pd
import os
import sys
import multiprocessing as mp


ds = sys.argv[1]
ds_dir = sys.argv[2]
job_id = sys.argv[3]
excel_file = '/scratch/tsour.s/Supplementary_Data_7.SAAP_coordinates.xlsx' # Supplementary Data Table 7 contains genomic locations of amino acid substitutions


# Load coordinates and convert to a set of tuples for O(1) lookup
def get_target_coords(excel_path, ds_filter):
    coord_data = pd.read_excel(excel_path)
    # Filter for the specific dataset
    ds_coord_data = coord_data[coord_data['Datasets'].str.contains(ds_filter, na=False)]
    # Store as a set of (chr, pos) strings for fast comparison
    return set(zip(ds_coord_data['chr'].astype(str), ds_coord_data['coor'].astype(str)))

# function to get the read depth of a list of gemonic coordinates
def extract_saap_reads_fast(id, ds, ds_dir, target_coords):
    input_path = os.path.join(ds_dir, id, f"{id}.per_base_coverage.txt")
    output_path = os.path.join(ds_dir, f"{id}_SAAP_base_coverage.txt")
    if not os.path.exists(input_path) or os.path.getsize(input_path) == 0:
        print(f"Input file missing or empty: {input_path}")
        return
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        print(f"Output already exists: {output_path}")
        return
    print(f"Processing {id}...")
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        outfile.write('CHROM\tPOS\tN reads\n')
        # Stream line by line (Memory efficient)
        for line in infile:
            parts = line.split('\t')
            # Assuming format: CHROM \t POS \t READS ...
            if (parts[0], parts[1]) in target_coords:
                outfile.write(line)


targets = get_target_coords(excel_file, ds)
extract_saap_reads_fast(job_id, ds, ds_dir, targets)
