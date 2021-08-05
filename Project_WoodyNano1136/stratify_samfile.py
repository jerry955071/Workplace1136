import pandas as pd
import sys
from WoodyNano import samtools

path_woodynano = sys.argv[1]
path_pychopper = sys.argv[2]
path_align_summary = sys.argv[3]

df = pd.read_csv(path_align_summary)
c = (df['read_length_woodynano'] - df['read_length_pychopper'])
cond1 = (c >= min(c)) & (c < 0)  # D
cond2 = c == 0  # C
cond3 = (c > 0) & (c < -min(c))  # B
cond4 = c >= -min(c)  # A


def main(
    read_names,
    path_out,
    path_sam={'woodynano':'', 'pychopper':''}
    ):

    def last_name(alignment):
        return alignment.qname.split('|')[-1]

    # read samfiles
    sam_woodynano = samtools.SAM.Import(
        fname=path_sam['woodynano'],
        get_primary=True
    )
    sam_pychopper = samtools.SAM.Import(
        fname=path_sam['pychopper'],
        get_primary=True
    )

    for i in sam_woodynano.alignment:
        if not last_name(i) in read_names:
            sam_woodynano.alignment.remove(i)

    for i in sam_pychopper.alignment:
        if not last_name(i) in read_names:
            sam_pychopper.alignment.remove(i)

    sam_woodynano.export(fname=path_out['woodynano'])
    sam_pychopper.export(fname=path_out['pychopper'])

    return

conds = [cond1, cond2, cond3, cond4]
groups = ['D', 'C', 'B', 'A']

for cond, group in zip(conds, groups):
    prefix_woodynano = path_woodynano.split('.sam')[0]
    prefix_pychopper = path_pychopper.split('.sam')[0]
    main(
        read_names=df[cond]['readname'],
        path_out={
            'woodynano':f"{prefix_woodynano}_group_{group}.sam",
            'pychopper':f"{prefix_pychopper}_group_{group}.sam"
            },
        path_sam={'woodynano':path_woodynano,'pychopper':path_pychopper}
        )

# end of script
print('Job end...')
