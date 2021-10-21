# %%
import pandas as pd
import numpy as np
import sys
#%%
# csv_path = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Stats/Summarized.csv'
csv_path = sys.argv[1]
df = pd.read_csv(csv_path, index_col=0)

a = np.array(df[df['software'] == 'woodynano']['read_length'])
b = np.array(df[df['software'] == 'pychopper']['read_length'])
c = a - b
# %%

cond1 = (c >= min(c)) & (c < 0)

tmp1 = df[df['software'] == 'woodynano'][cond1]
tmp2 = df[df['software'] == 'pychopper'][cond1]

pychopper_residual = np.array(tmp2['read_length']) - np.array(tmp1['read_length'])
pychopper_unmapped_residual = np.array(tmp2['softclipped_nucleotides']) - np.array(tmp1['softclipped_nucleotides'])

mapped_rate = 1 - sum(pychopper_unmapped_residual) / sum(pychopper_residual)

print(f'total residual = {sum(pychopper_residual)}')
print(f'unmapped_residual = {sum(pychopper_unmapped_residual)}')
print(f'mapped rate = {mapped_rate}')
# %%
