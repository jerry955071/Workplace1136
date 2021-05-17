#%%
from WoodyNano import samtools
import os
import sys


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
WoodyNano:{len(sam_woodynano.alignments)}\n\
Pychopper:{len(sam_pychopper.alignments)}'
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
    align_woodynano = dict_woodynano[name]
    align_pychopper = dict_pychopper[name]


# %%
