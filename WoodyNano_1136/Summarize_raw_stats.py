import pandas as pd
import os

# get csv filelist
csv_path = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Stats'
filelist = os.listdir(csv_path)
for i in filelist.copy():
    if not 'raw' in i:
        filelist.remove(i)

df = None
for name in filelist:
    tmp = pd.read_csv(f'{csv_path}/{name}')
    tmp['bioname'] = [name[:8]]*tmp.shape[0]
    df = pd.concat([df,tmp])

df = df.drop(['Unnamed: 0'], axis=1)
df.to_csv(f'{csv_path}/Summarized.csv')