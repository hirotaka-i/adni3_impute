# python code/make_bed_from_variant_list.py --variant_list /path/to/variant_list.tsv
import pandas as pd
import argparse
# Parse arguments
parser = argparse.ArgumentParser(description="Generate bed file for plink2 to slice the pfile")
parser.add_argument('--variant_list', help='Input variant list file path', required=True)   
args = parser.parse_args()
variant_list = args.variant_list

d=pd.read_csv(variant_list, sep='\t', names=['Chr', 'Start', 'End', 'Ref', 'Alt', 'rsid'])  
d_uniq=d[['Chr', 'Start', 'End']].drop_duplicates()

# Save the bed file
file_path = 'variants.bed1'
d_uniq.to_csv(file_path, index=False, header=False, sep='\t')
print(f'Bed file saved at {file_path}')