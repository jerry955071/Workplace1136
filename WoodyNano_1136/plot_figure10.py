#%%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


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
                              hue='Species_bio',
                              legend=False,
                              palette={'Ptr_bio1': 'darkorange',
                                       'Ptr_bio2': 'orange',
                                       'Egr_bio1': 'blue',
                                       'Egr_bio2': 'dodgerblue',
                                       'Lch_bio1': 'green',
                                       'Lch_bio2': 'limegreen'},
                              s=2,
                              ax=ax1)

    snsplt1.set(title=title,
                xlabel='Pychopper',
                ylabel='WoodyNano',
                xlim=lims,
                ylim=lims
                )

    snsplt2 = sns.histplot(x=data['Woodynano']-data['Pychopper'],
                           hue=data['Species_bio'],
                           ax=ax2,
                           palette={'Ptr_bio1': 'darkorange',
                                    'Ptr_bio2': 'orange',
                                    'Egr_bio1': 'blue',
                                    'Egr_bio2': 'dodgerblue',
                                    'Lch_bio1': 'green',
                                    'Lch_bio2': 'limegreen'},
                           fill=False,
                           legend=False,
                           element='poly',
                           binwidth=10)
    snsplt2.set(xlabel='WoodyNano-Pychopper', xlim=[-276, 276])

    fig.savefig(path)

    return

path_summarized_csv = '/Users/zhujiachen/Desktop/WoodyNano_Revision/Stats/Summarized.csv'
tmp = pd.read_csv(path_summarized_csv)

reshaped_df = {
    'Woodynano': [],
    'Pychopper': [],
    'Species_bio': []
}

reshaped_df['Woodynano'] = list(
    tmp[tmp['software'] == 'woodynano']['read_length'])
reshaped_df['Pychopper'] = list(
    tmp[tmp['software'] == 'pychopper']['read_length'])
reshaped_df['Species_bio'] = list(
    tmp[tmp['software'] == 'pychopper']['bioname'])

reshaped_df = pd.DataFrame(reshaped_df)

Figure10_modified(
    data=reshaped_df,
    path='/Users/zhujiachen/Desktop/WoodyNano_Revision/Figure10_revision.png',
    title=''
)
#%%