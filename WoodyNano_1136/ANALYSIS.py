#%%
from WoodyNano import samtools
import os
import sys
import pandas as pd

def last_name(alignment):
    return alignment.qname.split('|')[-1]

class logger(object):
    def __init__(self, path=None):
        super().__init__()
        self._path = path

    def new(self, times=None):
        with open(self._path, 'w'):
            os.utime(self._path, times)
        return

    def append(self, lines=None):
        with open(self._path, 'a') as o:
            o.write(lines+'\n')
        return

# path_woodynano = sys.argv[1]
# path_pychopper = sys.argv[2]
# path_logfile = sys.argv[3]
path_woodynano = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_woodynano.sam'
path_pychopper = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_full_length_pychopper.sam'
path_logfile = '/Users/zhujiachen/Desktop/remap/Ptr_bio1_sam_compare.log'
out_table = '/Users/zhujiachen/Desktop/remap/Stats/Ptr_bio1_stats.csv'

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
        both_mapped.pop(name)

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


# %%
