#%%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
#%%
csv_path = sys.argv[1]
figure_path = sys.argv[2]

def Figure10_modified(data, path, title, lims=[0, 6000]):
    sns.set_theme(
        context='paper',
        font_scale=1.2,
        font='sans serif'
    )
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 12))
    ax1.plot(lims, lims, ls=':', lw=2)
    ax2.axvline(ls=':', color='black', lw=2)
    snsplt1 = sns.scatterplot(data=data,
                              x='Pychopper',
                              y='Woodynano',
                              legend=False,
                              s=2,
                              ax=ax1)

    snsplt1.set(title=title,
                xlabel='Pychopper',
                ylabel='WoodyNano',
                xlim=lims,
                ylim=lims
                )

    snsplt2 = sns.histplot(x=data['Woodynano']-data['Pychopper'],
                           ax=ax2,
                           fill=False,
                           legend=False,
                           element='poly',
                           binwidth=10)
    snsplt2.set(xlabel='WoodyNano-Pychopper', xlim=[-276, 276])

    fig.savefig(path)

    return


path_summarized_csv = csv_path
tmp = pd.read_csv(path_summarized_csv)

df = tmp[['read_length_woodynano', 'read_length_pychopper']]
df.columns = ['Woodynano', 'Pychopper']

#%%
Figure10_modified(
    data=df,
    path=figure_path,
    title=''
)
#%%
