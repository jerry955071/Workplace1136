# %%
import sys
from WoodyNano import samtools
import os
import pandas as pd
import numpy as np

def last_name(alignment):
    return alignment.qname.split('|')[-1]

path_woodynano = sys.argv[1]
path_pychopper = sys.argv[2]
path_align_summary = sys.argv[3]

# path_woodynano = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/Ptr_bio1_full_length_woodynano.sam'
# path_pychopper = '/Users/zhujiachen/Desktop/WoodyNano_Revision/remap/Ptr_bio1_full_length_pychopper.sam'

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
df.to_csv(path_align_summary, index=False)
