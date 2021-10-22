import pandas as pd
import numpy as np
import sys

# Usage:
# python [scirpt] [2in1.csv] > [output]

# [2in1.csv] file header:
# readname,read_length_woodynano,read_length_pychopper,
# mapped_read_length_woodynano,mapped_read_length_pychopper,
# matching_nucleotides_woodynano,matching_nucleotides_pychopper,
# softclipped_nucleotides_woodynano,softclipped_nucleotides_pychopper,
# clipping_left_woodynano,clipping_left_pychopper,
# clipping_right_woodynano,clipping_right_pychopper,
# ref_start_woodynano,ref_start_pychopper,
# ref_end_woodynano,ref_end_pychopper,
# numbers_of_exon_woodynano,numbers_of_exon_pychopper

csv_path = sys.argv[1]
df = pd.read_csv(csv_path, index_col=0)

a = np.array(df["read_length_woodynano"])
b = np.array(df["read_length_pychopper"])
c = a - b

cond1 = (c >= min(c)) & (c < 0)

# tmp1 = df["read_length_woodynano"][cond1]
# tmp2 = df["read_length_pychopper"][cond1]

pychopper_residual = np.array(df["read_length_pychopper"][cond1]) - np.array(df["read_length_woodynano"][cond1])
pychopper_unmapped_residual = np.array(df["softclipped_nucleotides_pychopper"][cond1]) - np.array(df['softclipped_nucleotides_woodynano'][cond1])

mapped_rate = 1 - sum(pychopper_unmapped_residual) / sum(pychopper_residual)

print(f'total residual = {sum(pychopper_residual)}')
print(f'unmapped_residual = {sum(pychopper_unmapped_residual)}')
print(f'mapped rate = {mapped_rate}')
# %%
