import pandas as pd
import os
import sys

# get csv filelist
# csv_path = '/home/woodformation/Processing_data/CCC/Processed_data_WoodyNano/Stats'
try:
    csv_path = sys.argv[1]
except:
    raise Exception('Path of stats were required')
    
filelist = os.listdir(csv_path)
for i in filelist.copy():
    if not 'all' in i:
        filelist.remove(i)

df = None
for name in filelist:
    tmp = pd.read_csv(f'{csv_path}/{name}')
    tmp['bioname'] = [name[:8]]*tmp.shape[0]
    df = pd.concat([df,tmp])

df.to_csv(f'{csv_path}/Summarized.csv', index=False)
