#%%
from WoodyNano import seqtools
import sys
import numpy as np

def last_name(readname):
    return readname.split('|')[-1].split(' ')[0]

args = sys.argv

# raw_fastq = seqtools.SeqFastq.Import(fname=args[1])
# woodynano_fastq = seqtools.SeqFastq.Import(fname=args[2])
# pychopper_fastq = seqtools.SeqFastq.Import(fname=args[3])

raw_fastq = seqtools.SeqFastq.Import(
    fname='/home/woodformation/Processing_data/CCC/Rawdata_WoodyNano/Fastq/Raw/Kit109/Egr_bio1_kit109.fastq')
woodynano_fastq = seqtools.SeqFastq.Import(
    fname='/home/woodformation/Processing_data/CCC/Rawdata_WoodyNano/Fastq/Trimmed/Woodynano/Egr_bio1_full_length_woodynano.fastq')
pychopper_fastq = seqtools.SeqFastq.Import(
    fname='/home/woodformation/Processing_data/CCC/Rawdata_WoodyNano/Fastq/Trimmed/Pychopper/Egr_bio1_full_length_pychopper.fastq')

for i in list(woodynano_fastq.keys()):
    woodynano_fastq[last_name(i)] = woodynano_fastq.pop(i)

for i in list(pychopper_fastq.keys()):
    pychopper_fastq[last_name(i)] = pychopper_fastq.pop(i)

woodynano = woodynano_fastq.keys()
pychopper = pychopper_fastq.keys()


intersect = set(woodynano) & set(pychopper)
woodynano_uni = set(woodynano) - intersect
pychopper_uni = set(pychopper) - intersect

#rownames = ['intersect','woodynano unique','pychopper unique']
#for name_set, rowname in zip([intersect, woodynano_uni, pychopper_uni], rownames):
#    avg_qscores = list()
#    for name in name_set:
#        avg_qscores.append(raw_fastq[name].mean_q())
#
#    print('%s: %.2f %.2f %.2f' % (rowname, np.quantile(avg_qscores, .25), \
#            np.mean(avg_qscores), np.quantile(avg_qscores, .75)))
#


#%%
for name_set, rowname in zip([intersect, pychopper_uni],['pychopper_intersect','pychopper_unique']):
    avg_qscores = list()
    for name in name_set:
        avg_qscores.append(pychopper_fastq[name].mean_q())

    print('%s: %.2f' % (rowname, np.mean(avg_qscores)))

#%%

for name_set, rowname in zip([intersect, woodynano_uni], ['woodynano_intersect', 'woodynano_unique']):
    avg_qscores = list()
    for name in name_set:
        avg_qscores.append(woodynano_fastq[name].mean_q())

    print('%s: %.2f' % (rowname, np.mean(avg_qscores)))


# %%
