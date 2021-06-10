#%%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
#%%
csv_path = sys.argv[1]
figure10_path = sys.argv[2]

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


path_summarized_csv = csv_path
tmp = pd.read_csv(path_summarized_csv)

reshaped_df = tmp[['bioname', 'read_length_woodynano', 'read_length_pychopper']]
slice_1 = reshaped_df[reshaped_df['bioname'] == 'Egr_bio1']
slice_2 = reshaped_df[reshaped_df['bioname'] == 'Egr_bio2']
slice_3 = reshaped_df[reshaped_df['bioname'] == 'Ptr_bio1']
slice_4 = reshaped_df[reshaped_df['bioname'] == 'Ptr_bio2']
slice_5 = reshaped_df[reshaped_df['bioname'] == 'Lch_bio1']
slice_6 = reshaped_df[reshaped_df['bioname'] == 'Lch_bio2']

sorted_df = pd.concat([slice_1, slice_2, slice_3, slice_4, slice_5, slice_6])
sorted_df.columns = ['Species_bio', 'Woodynano', 'Pychopper']

#%%
Figure10_modified(
    data=sorted_df,
    path=figure10_path,
    title=''
)
#%%
