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
cond_bio = df['bioname'] == path_woodynano.split('/')[-1].split('_full')[0]

bio_name = path_woodynano.split('/')[-1].split('_full')[0]
print(bio_name)

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

    new_woodynano = samtools.SAM()
    new_woodynano.set_header(sam_woodynano.get_header())
    new_pychopper = samtools.SAM()
    new_pychopper.set_header(sam_pychopper.get_header())

    for i in sam_woodynano.alignment:
        if last_name(i) in read_names:
            new_woodynano.alignment.append(i)

    for i in sam_pychopper.alignment:
        if last_name(i) in read_names:
            new_pychopper.alignment.append(i)

    print(len(new_woodynano.alignment), len(new_pychopper.alignment))

    new_woodynano.export(fname=path_out['woodynano'])
    new_pychopper.export(fname=path_out['pychopper'])

    return

conds = [cond1, cond2, cond3, cond4]
groups = ['D', 'C', 'B', 'A']

for cond, group in zip(conds, groups):
    print(group, sum(cond))
    prefix_woodynano = path_woodynano.split('.sam')[0]
    prefix_pychopper = path_pychopper.split('.sam')[0]

    print(f"{prefix_woodynano}_group_{group}.sam", f"{prefix_pychopper}_group_{group}.sam")
    
    read_names = set(df[cond & cond_bio]['readname'])

    main(
        read_names,
        path_out={
            'woodynano':f"{prefix_woodynano}_group_{group}.sam",
            'pychopper':f"{prefix_pychopper}_group_{group}.sam"
            },
        path_sam={'woodynano':path_woodynano,'pychopper':path_pychopper}
        )

# end of script
print('Job end...')
