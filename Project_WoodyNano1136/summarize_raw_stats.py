import pandas as pd
import os
import sys

# Usage:
# python [script] [fout] [filenames] [samplenames]

# get csv filelist
# csv_path = '/home/woodformation/Processing_data/CCC/Processed_data_WoodyNano/Stats'
# try:
#     csv_path = sys.argv[1]
# except:
#     raise Exception('Path of stats were required')
    
# filelist = os.listdir(csv_path)
# for i in filelist.copy():
#     if not 'all' in i:
#         filelist.remove(i)

# df = None
# for name in filelist:
#     tmp = pd.read_csv(f'{csv_path}/{name}')
#     tmp['bioname'] = [name[:8]]*tmp.shape[0]
#     df = pd.concat([df,tmp])

# df.to_csv(f'{csv_path}/Summarized.csv', index=False)

args = sys.argv[1:]
fout = args.pop(0)
length = len(args)
filenames = args[0:length/2]
samplenames = args[length/2:]
print(filenames)
print(samplenames)

df_list = [pd.read_csv(i) for i in filenames]
for df, name in zip(df_list, samplenames):
    df["sample"] = [name] * df.shape[0]

df_ttl = pd.concat(df_list)
df_ttl.to_csv(fout)
    
