# usage:
# python code/adjust_score_file.py\
#  --pvar pvar\
#  --rsid_ref rsid_ref\
#  --score_file score_file\
#  --out out

import pandas as pd
import numpy as np
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Adjust PRS score file")
parser.add_argument('--pvar', help='Input pvar file path', required=True)
parser.add_argument('--rsid_ref', help='Input rsid reference file path', required=True)
parser.add_argument('--score_file', help='Input PRS file path', required=True)
parser.add_argument('--out', help='Output file path', required=True)
args = parser.parse_args()
pvar_path = args.pvar
rsid_path = args.rsid_ref
score_file = args.score_file
out = args.out

def find_header_row(filename):
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith('#CHROM'):
                return i
    return None

# read pvar
header_row = find_header_row(pvar_path)
pvar = pd.read_csv(pvar_path, sep='\s+', header=header_row)

# Add rsid to the pvar
d_rsID_uniq=pd.read_csv(rsid_path, sep='\t', names=['Chr', 'Start', 'End', 'Ref', 'Alt', 'rsid']).drop_duplicates(subset=['Chr', 'Start', 'rsid'])
r = pd.merge(pvar, d_rsID_uniq, left_on=['#CHROM', 'POS'], right_on=['Chr', 'Start'], how='left')


# adjust score file
print(f'Start adjusting for {score_file}')
t = pd.read_csv(score_file, sep='\t', header=None, names=['rsid', 'A1', 'Beta'])
print(f'\nINPUT score file rows: {t.shape[0]}')
d=pd.merge(r, t, on=['rsid'], how='inner')
print(f'Missing\n{np.setdiff1d(t.rsid, d.rsid)}')

# duplicates
d_dups=d[d.duplicated(subset=['rsid'], keep=False)].sort_values(by='rsid')
if d_dups.shape[0]>0:
    print(f'score file merged into the (split) multiple alleles: need to resovle {d_dups.shape[0]} duplicated entries')
    print('first remmedy: A1 is ALT --> taking ID with corresponding ALT allele')
    d_dups_1=d_dups[d_dups.A1==d_dups.ALT]
    print(f'first remedy applied: {d_dups_1.shape[0]}')
    print('second remmedy: A1 is REF --> take the one with more freq entries')
    d_dups_2=d_dups[d_dups.A1==d_dups.REF]

    # d_dups_2=d_dups[d_dups.A1==d_dups.REF].sort_values('ALT_FREQS').drop_duplicates(subset=['rsid'], keep='last')
    print(f'second remedy applied: {d_dups_2.shape[0]}')

    # concat 1 and 2
    d_dup_resolved=pd.concat([d_dups_1, d_dups_2], axis=0, ignore_index=True)

    # non-duplicated entries
    d_nodups=d[~d.duplicated(subset=['rsid'], keep=False)]

    # combine together
    d=pd.concat([d_nodups, d_dup_resolved], axis=0, ignore_index=True).sort_values(['#CHROM', 'POS'])

# allele adjust
if 'AD_Bellenguez' in score_file:
    # adjust the score file
    if 'chr4:993555' in d.ID.values:
        if d[d.ID=='chr4:993555']['REF'].values[0]=='GAGTT':
            d.loc[d.ID=='chr4:993555', 'A1']='TAGTT'
            print('AD_Bellenguez: chr4:11023507 T is adjusted to TAGTT because the REF is GAGTT')

print(f'OUTPUT: adjusted score file rows: {d.shape[0]}')
# save the score file
d[['ID', 'A1', 'Beta']].to_csv(out, index=False, header=False, sep='\t')