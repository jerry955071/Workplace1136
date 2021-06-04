from WoodyNano import seqtools
from WoodyNano import logger
import sys


def last_name(read_info):
    tmp = read_info.split('|')[-1]
    return tmp.split(' ')[0]

fastq_woodynano = sys.argv[1]
fastq_pychopper = sys.argv[2]
output_file = sys.argv[3]

readname_woodynano = list()
readname_pychopper = list()

with open(fastq_woodynano) as f:
    read = seqtools.SeqFastq.read(f)
    while read:
        readname_woodynano.append(last_name(read.info))
        read = seqtools.SeqFastq.read(f)

with open(fastq_pychopper) as f:
    read = seqtools.SeqFastq.read(f)
    while read:
        readname_pychopper.append(last_name(read.info))
        read = seqtools.SeqFastq.read(f)

intersect = set(readname_woodynano) & set(readname_pychopper)
len_inter = len(intersect)
w_unique = len(readname_woodynano) - len_inter
q_unique = len(readname_pychopper) - len_inter

out_file = logger.logger(output_file)
out_file.new()
out_file.append('w_unique\tintersect\tp_unique')
out_file.append(f'{w_unique}\t{len_inter}\t{q_unique}')

