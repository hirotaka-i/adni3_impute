import pandas as pd
import os
data_dir = os.getenv('DATA_DIR')
temp_dir = os.getenv('TEMP_DIR')

variant_file = f'{data_dir}/PRS/rsid_of_interest_chrpos.txt'
cmd = f"python ../make_bed_from_variant_list.py --variant_list {variant_file}"
