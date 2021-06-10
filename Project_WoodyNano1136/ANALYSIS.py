#%%
from WoodyNano import samtools
from WoodyNano import logger
import os
import sys
import pandas as pd
import numpy as np

def last_name(alignment):
    return alignment.qname.split('|')[-1]

# path_woodynano = sys.argv[1]
# path_pychopper = sys.argv[2]
# path_logfile = sys.argv[3]
path_woodynano = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_woodynano.sam'
path_pychopper = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_pychopper.sam'
path_logfile = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_sam_compare.log'
out_table = '/Users/zhujiachen/Desktop/remap/Stats/Ptr_bio1_stats.csv'
out_table_summarized = '/Users/zhujiachen/Desktop/remap/Stats/Ptr_bio1_stats_summarized.csv'


logfile = logger(path=path_logfile)
logfile.new()
logfile.append(f'Compare the alignment results between:\n{path_woodynano} and\n{path_pychopper}')

sam_woodynano = samtools.SAM.Import(
    fname=path_woodynano,
    get_primary=True
)
sam_pychopper = samtools.SAM.Import(
    fname=path_pychopper,
    get_primary=True
)

logfile.append(
    f'Numbers of primary alignments:\n\
WoodyNano:{len(sam_woodynano.alignment)}\n\
Pychopper:{len(sam_pychopper.alignment)}'
)

dict_woodynano = {last_name(i):i for i in sam_woodynano.alignment}
dict_pychopper = {last_name(i):i for i in sam_pychopper.alignment}

both_mapped = set(dict_woodynano.keys()) & set(dict_pychopper.keys())
logfile.append(f'Numbers of both mapped primary alignments:{both_mapped.__len__()}')

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

logfile.append(f'Numbers of both mapped/same loci/primary alignments:{same_loci.__len__()}')


stats = {
    'readname':[],
    'software':[],
    'read_length':[],
    'mapped_read_length':[],
    'matching_nucleotides':[],
    'softclipped_nucleotides':[]
}

for name in same_loci:
    for software, dictionary in zip(['woodynano','pychopper'],[dict_woodynano, dict_pychopper]):
        stats['readname'].append(name)
        stats['software'].append(software)
        stats['read_length'].append(dictionary[name].read_length)
        stats['mapped_read_length'].append(dictionary[name].mapped_read_length)
        stats['matching_nucleotides'].append(dictionary[name].matching_nucleotides)
        stats['softclipped_nucleotides'].append(sum(dictionary[name].clipping))

df = pd.DataFrame(stats)
df.to_csv(out_table)
logfile.append(f'Stats output (raw) at:{out_table}')

# %%
table1 = {
    'parameters':[],
    'software':[],
    'mean':[],
    'quantile_0.25':[],
    'median':[],
    'quantile_0.75':[]
}


for parameters in 'read_length','mapped_read_length','matching_nucleotides','softclipped_nucleotides':
    for software in 'woodynano','pychopper':
        table1['software'].append(software)
        table1['parameters'].append(parameters)
        table1['mean'].append(
            np.mean(
                df[ df['software'] == software ][parameters]
                )
            )
        quantiles = np.quantile(
                a=df[ df['software'] == software ][parameters],
                q=[0.25, 0.5, 0.75]
                )
        table1['quantile_0.25'].append(quantiles[0])
        table1['median'].append(quantiles[1])
        table1['quantile_0.75'].append(quantiles[2])

for col in 'mean','quantile_0.25', 'median', 'quantile_0.75':
    table1[col] = [round(i,1) for i in table1[col]]

logfile.append(f'Stats output (summarized) at:{out_table_summarized}')
pd.DataFrame(table1).to_csv(out_table_summarized)

# %%
def try_atom():
    print('This line is edit with atoms')
