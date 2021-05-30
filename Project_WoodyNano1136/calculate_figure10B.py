# %%
import pandas as pd
import numpy as np
#%%
csv_path = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Stats/Summarized.csv'
df = pd.read_csv(csv_path, index_col=0)


# %%
a = np.array(df[df['software'] == 'woodynano']['read_length'])
b = np.array(df[df['software'] == 'pychopper']['read_length'])
c = a - b
# %%

min(c), 0, -min(c)

cond1 = (c >= min(c)) & (c < 0)
cond2 = c == 0
cond3 = (c > 0) & (c < -min(c))
cond4 = c >= -min(c)

print(*[f'{sum(cond)/len(c)*100}\n' for cond in [cond1, cond2, cond3, cond4]])
# %%

# min(c), 0, -min(c)...
# 96.74998756186145
#  0.26017659230123125
#  1.9312770895341114
#  1.0585587563032097