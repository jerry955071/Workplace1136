# %%
import sys
from WoodyNano import samtools
import os
import pandas as pd
import numpy as np

def last_name(alignment):
    return alignment.qname.split('|')[-1]

# path_woodynano = sys.argv[1]
# path_pychopper = sys.argv[2]

path_woodynano = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/Ptr_bio1_full_length_woodynano.sam'
path_pychopper = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/Ptr_bio1_full_length_pychopper.sam'

sam_woodynano = samtools.SAM.Import(
    fname=path_woodynano,
    get_primary=True
)
sam_pychopper = samtools.SAM.Import(
    fname=path_pychopper,
    get_primary=True
)


dict_woodynano = {}
dict_pychopper = {}

for i in sam_woodynano.alignment:
    if not last_name(i) in dict_woodynano.keys():
        dict_woodynano[last_name(i)] = i
for i in sam_pychopper.alignment:
    if not last_name(i) in dict_pychopper.keys():
        dict_pychopper[last_name(i)] = i


both_mapped = set(dict_woodynano.keys()) & set(dict_pychopper.keys())

for i in set(dict_woodynano.keys()) - both_mapped:
    dict_woodynano.pop(i)

for i in set(dict_pychopper.keys()) - both_mapped:
    dict_pychopper.pop(i)

# Check if on the same loci
for name in list(both_mapped):
    dict_woodynano[name].cal_all_stats(keep=True)
    dict_pychopper[name].cal_all_stats(keep=True)


for name in list(both_mapped):
    ref_span_woodynano = dict_woodynano[name].reference_span
    ref_span_pychopper = dict_pychopper[name].reference_span
    if not min(ref_span_woodynano[-1],ref_span_pychopper[-1]) > max(ref_span_woodynano[0],ref_span_pychopper[0]):
        both_mapped.remove(name)

same_loci = both_mapped

#%%

stats = {
    'readname':[],
    'read_length_woodynano':[],
    'read_length_pychopper':[],
    'mapped_read_length_woodynano':[],
    'mapped_read_length_pychopper':[],
    'matching_nucleotides_woodynano':[],
    'matching_nucleotides_pychopper':[],
    'softclipped_nucleotides_woodynano':[],
    'softclipped_nucleotides_pychopper':[],
    'clipping_left_woodynano':[],
    'clipping_left_pychopper':[],
    'clipping_right_woodynano':[],
    'clipping_right_pychopper':[],
    'ref_start_woodynano':[],
    'ref_start_pychopper':[],
    'ref_end_woodynano':[],
    'ref_end_pychopper':[],
    'numbers_of_exon_woodynano':[],
    'numbers_of_exon_pychopper':[]
}

for name in same_loci:
    stats['readname'].append(name)
    for software, dictionary in zip(['woodynano','pychopper'],[dict_woodynano, dict_pychopper]):
        num_extron = dictionary[name].cigar_summary['extron'].__len__()
        num_exon = (num_extron - 1) / 2 + 1
        stats[f'read_length_{software}'].append(dictionary[name].read_length)
        stats[f'mapped_read_length_{software}'].append(dictionary[name].mapped_read_length)
        stats[f'matching_nucleotides_{software}'].append(dictionary[name].matching_nucleotides)
        stats[f'softclipped_nucleotides_{software}'].append(sum(dictionary[name].clipping))
        stats[f'clipping_left_{software}'].append(dictionary[name].clipping[0])
        stats[f'clipping_right_{software}'].append(dictionary[name].clipping[1])
        stats[f'ref_start_{software}'].append(dictionary[name].reference_span[0])
        stats[f'ref_end_{software}'].append(dictionary[name].reference_span[1])
        stats[f'numbers_of_exon_{software}'].append(int(num_exon))

df = pd.DataFrame(stats)
df.to_csv('/Users/zhujiachen/Desktop/WoodyNano_Revision/alignment_summary.csv', index=False)
# %%
# stratified 'df' with the difference in read length 
df_w = df[df['software'] == 'woodynano']
df_p = df[df['software'] == 'pychopper']

a = np.array(df_w['read_length'])
b = np.array(df_p['read_length'])
c = a - b

cond_d = (c >= min(c)) & (c < 0)
cond_c = c == 0
cond_b = (c > 0) & (c < -min(c))
cond_a = c >= -min(c)

def summary(array):
    return f'mean:\t{np.mean(array)}\
        \nquantile 0:\t{np.quantile(array, 0)}\
        \nquantile 0.25:\t{np.quantile(array, .25)}\
        \nquantile 0.5:\t{np.quantile(array, .5)}\
        \nquantile 0.75:\t{np.quantile(array, .75)}\
        \nquantile 1.0:\t{np.quantile(array, 1.0)}\
        '

def count_numbers(array):
    numbers = set(array)
    out_str = str()
    for n in numbers:
        out_str += f'WoodyNano identified {n} more exon(s):\t{sum(array == n)}\n'
    
    return out_str


conds = cond_a, cond_b, cond_c, cond_d
descriptions = ['group A', 'group B', 'group C', 'group D']
print('Numbers of exon (W-P)')
for cond,description in zip(conds,descriptions):
    nexon_w = np.array(df_w[cond]['numbers_of_exon'])
    nexon_p = np.array(df_p[cond]['numbers_of_exon'])
    diff_nexon = nexon_w - nexon_p
    # print(description,summary(diff_nexon),'', sep='\n')
    print(f'{description}, {sum(cond)} reads', count_numbers(diff_nexon), sep='\n')

# %%
nexon_w = np.array(df_w[cond_a]['numbers_of_exon'])
nexon_p = np.array(df_p[cond_a]['numbers_of_exon'])
diff_nexon = nexon_w - nexon_p
groupA_diff_nexon_woodynano_plus1 = df_w[cond_a][diff_nexon == 1]['readname']

new_sam_woodynano = samtools.SAM()
new_sam_pychopper = samtools.SAM()

new_sam_woodynano.set_header(sam_woodynano.get_header())
new_sam_pychopper.set_header(sam_woodynano.get_header())

new_sam_woodynano.alignment = list()
new_sam_pychopper.alignment = list()

for i in groupA_diff_nexon_woodynano_plus1:
    new_sam_woodynano.alignment.append(dict_woodynano[i])
    new_sam_pychopper.alignment.append(dict_pychopper[i])

new_sam_woodynano.export(f'/Users/zhujiachen/Desktop/WoodyNano_Revision/Group_SAM/groupA_diff_nexon_woodynano_plus1_woodynano.sam')
new_sam_pychopper.export(f'/Users/zhujiachen/Desktop/WoodyNano_Revision/Group_SAM/groupA_diff_nexon_woodynano_plus1_pychopper.sam')
# %%
