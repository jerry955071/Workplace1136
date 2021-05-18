import seaborn as sns
import pandas as pd
import os

# read csv
csv_path = str()
filelist = os.listdir(csv_path)

for i in filelist.copy():
    if 'summarized' in i:
        filelist.pop(i)

for clean in filelist:
    bioname = clean.split('/')[-1][:7]
    df = pd.read_csv(csv_path)
    nrow = df.shape[0]
    df['bioname'] = [bioname]*nrow

    try:
        pd.concat([tmp, df])
    except NameError:
        tmp = df
    