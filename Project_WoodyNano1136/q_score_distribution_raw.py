from WoodyNano import seqtools
import sys
import numpy as np

def last_name(readname):
    return readname.split('|')[-1].split(' ')[0]

args = sys.argv

raw_fastq = seqtools.SeqFastq.Import(fname=args[1])
woodynano = [last_name(i) for i in seqtools.SeqFastq.Import(fname=args[2]).keys()]
pychopper = [last_name(i) for i in seqtools.SeqFastq.Import(fname=args[3]).keys()]


intersect = set(woodynano) & set(pychopper)
woodynano_uni = set(woodynano) - intersect
pychopper_uni = set(pychopper) - intersect

rownames = ['intersect','woodynano unique','pychopper unique']
for name_set, rowname in zip([intersect, woodynano_uni, pychopper_uni], rownames):
    avg_qscores = list()
    for name in name_set:
        avg_qscores.append(raw_fastq[name].mean_q())

    print('%s: %.2f %.2f %.2f' % (rowname, np.quantile(avg_qscores, .25), \
            np.mean(avg_qscores), np.quantile(avg_qscores, .75)))




