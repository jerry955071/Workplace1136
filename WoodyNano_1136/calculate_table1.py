#%%
import pandas as pd
import numpy as np
#%%
csv_path = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Stats/Summarized.csv'
df = pd.read_csv(csv_path)
# %%
stats = ['read_length','mapped_read_length','matching_nucleotides','softclipped_nucleotides']
table1 = {
    'Bio':[],
    'Statistic': [],
    'Woodynano':[],
    'Pychopper':[]
}

for stat in stats:
    for bio in 'Egr_bio1','Egr_bio2','Ptr_bio1','Ptr_bio2','Lch_bio1','Lch_bio2':
        cond = df['bioname'] == bio
        cond_w = df['software'] == 'woodynano'
        cond_p = df['software'] == 'pychopper'
        table1['Bio'].append(bio)
        table1['Statistic'].append(stat)
        table1['Woodynano'].append(np.mean(df[cond & cond_w][stat]))
        table1['Pychopper'].append(np.mean(df[cond & cond_p][stat]))

# %%
pd.DataFrame(table1).to_csv('/Users/zhujiachen/Desktop/WoodyNano_Revision/Table1.csv',index=False)
# %%
