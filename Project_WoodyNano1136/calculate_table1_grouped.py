#%%
import pandas as pd
import numpy as np
import sys
#%%
csv_path = sys.argv[1]
path_table1 = sys.argv[2]

# csv_path = '/home/woodformation/Processing_data/CCC/Processed_data_WoodyNano/Stats/Summarized.csv'
df = pd.read_csv(csv_path)
# %%
c = df['read_length_woodynano'] - df['read_length_pychopper']
cond1 = (c >= min(c)) & (c < 0)  # D		96.6880962193937
cond2 = c == 0  # C						0.26052399618539923
cond3 = (c > 0) & (c < -min(c))  # B		1.9632558881130664
cond4 = c >= -min(c)  # A 				1.0881238963078395

for cond, group in zip([cond1, cond2, cond3, cond4], ['D', 'C', 'B', 'A']):

    stats = ['read_length','mapped_read_length','matching_nucleotides','softclipped_nucleotides']
    table1 = {
        'Bio':[],
        'Statistic': [],
        'Woodynano': [],
        'Pychopper': [],
        'Difference': [],
        'Average of difference': []
    }

    for stat in stats:
        diff_list = []
        diff_per_list = []
        for bio in 'Egr_bio1','Egr_bio2','Ptr_bio1','Ptr_bio2','Lch_bio1','Lch_bio2':
            cond0 = df['bioname'] == bio
            mean_w = np.mean(df[cond0 & cond][f'{stat}_woodynano'])
            mean_p = np.mean(df[cond0 & cond][f'{stat}_pychopper'])
            diff = mean_w - mean_p
            diff_per = diff/mean_w*100
            diff_list.append(diff)
            diff_per_list.append(diff_per)

            table1['Bio'].append(bio)
            table1['Statistic'].append(stat)
            table1['Woodynano'].append('%.1f' % mean_w)
            table1['Pychopper'].append('%.1f' % mean_p)
            table1['Difference'].append('%.1f (%.1f%%)' % (diff, diff_per))

        table1['Average of difference'] += ['%.1f (%.1f%%)' % (np.mean(diff_list), np.mean(diff_per_list)), '', '', '', '', '']

    # %%
    pd.DataFrame(table1).to_csv(f'{path_table1}/table1_group_{group}.csv', index=False)
    # %%
