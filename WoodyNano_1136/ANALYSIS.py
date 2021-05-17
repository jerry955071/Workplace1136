#%%
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from WoodyNano.Fastq import SeqFastq
from WoodyNano.paftools import PAF
'''
analysis procedure:
nanopore-raw.fastq --> full-length.fastq --> minimap-output.paf --> mapped-reads.paf
								|-> unmapped-reads.fasta --> blast-output.csv
'''


dirs = {}
dirs['full-length'] = '/home/woodformation/Processing_data/CCC/WoodyNano_raw/Fastq/Trimmed/'
dirs['paf'] = '/home/woodformation/Processing_data/CCC/WoodyNano_raw/PAF/'
dirs['blast'] = '/home/woodformation/Processing_data/CCC/WoodyNano_raw/BLAST/'
dirs['stats'] = '/home/woodformation/Processing_data/CCC/WoodyNano_raw/Stats/'


class logger(object):
    def __init__(self, path=None):
        super().__init__()
        self._path = path

    def new(self, times=None):
        with open(self._path, 'w'):
            os.utime(self._path, times)
        return

    def append(self, lines=None):
        with open(self._path, 'a') as o:
            o.write(lines+'\n')
        return


def get_readname(alignment):
    return alignment.query_name.split('|')[-1]


def get_path(parent, pattern):
    """
    parent: parent path
    pattern: pattern of file name
    search for file name matching the pattern in the parent directory, 
    return list[]: full path ("parent path/fname")
    """
    match = []
    for fname in os.listdir(path=parent):
        if re.search(pattern, fname):
            match.append(parent+fname)
    if len(match) == 1:
        match = str(match[0])
        return match
    pass


def statisics(software, species, bio, dirs):
    """
    Do not use
    """
    print('Do statistics')
    out_log = logger(path=f"{dirs['stats']}{software}/{species}_{bio}.stats")
    out_log.new()

    # read fastq
    # get file path to 'foo'
    foo = get_path(parent=f"{dirs['full-length']}{software}/",
                   pattern=f'{species}_{bio}'
                   )
    fq = SeqFastq.Import(fname=foo,
                         ftype='Fastq')
    # read paf
    # get file path
    foo = get_path(parent=f"{dirs['paf']}{software}/",
                   pattern=f'{species}_{bio}'
                   )
    paf = PAF.Import(fname=foo).get_primary_paf()

    # stats of full-length output
    strands = np.array([fq[n].strand for n in fq.keys()])
    out_log.append(lines='Full-length results:')
    out_log.append(lines=f'Numbers of full-length: {fq.__len__()}')
    out_log.append(
        lines=f'Total read length: {np.sum([fq[n].rlen for n in fq.keys()])}')
    out_log.append(lines='Strand(+/-): {rate:.4f}'.format(
        rate=np.sum(strands == '+')/np.sum(strands == '-')))

    # stats of mapping results
    mapped_strand = np.array([align.strand for align in paf.alignment])
    out_log.append(lines='\nMapping results:')
    out_log.append(lines='Mapping rate: {rate:.4f}'.format(
        rate=len(set([get_readname(i) for i in paf.alignment]))/fq.__len__())
    )
    out_log.append(
        lines=f'Mapped to forward strand: {np.mean(mapped_strand=="+")/len(paf.alignment)}')
    out_log.append(
        lines=f'Average effective read length: {np.mean([i.effective_read_length() for i in paf.alignment])}')
    out_log.append(
        lines=f'Average align block length: {np.mean([i.align_block_len for i in paf.alignment])}')
    out_log.append(
        lines=f'Average residual matched: {np.mean([i.align_block_len for i in paf.alignment])}')
    out_log.append(
        lines=f'Average clipped bases: {np.mean([i.clipped_base() for i in paf.alignment])}')
    #out_log.append(lines='Blasted rate: {rate:.4f}'.format(
    #    rate=len(set(blast['qaccver']))/fq.__len__())
    #)

    return


def compare(species, bio, dirs):
    print('Do compare')

    out_log = logger(
        path=f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/Stats/Comparison/{species}_{bio}.tsv')
    out_log.new()

    fq = {}
    paf = {}
    blast = {}
    for software in ['Woodynano', 'Pychopper']:
        # read fastq
        # get file path to 'foo'
        foo = get_path(parent=f"{dirs['full-length']}{software}/",
                       pattern=f'{species}_{bio}'
                       )
        fq[software] = SeqFastq.Import(fname=foo,
                                       ftype='Fastq')
        # read paf
        # get file path
        foo = get_path(parent=f"{dirs['paf']}{software}/",
                       pattern=f'{species}_{bio}'
                       )
        tmp = PAF.Import(fname=foo).get_primary_paf()
        paf[software] = {get_readname(i): i for i in tmp.alignment}
       
    # compare full-length output:
    out_log.append(lines='Compare full-length output')
    wd = set(fq['Woodynano'].keys())
    py = set(fq['Pychopper'].keys())
    inter = wd & py
    out_log.append('W.unique\tIntersection\tP.unique')
    out_log.append(f'{len(wd-inter)}\t{len(inter)}\t{len(py-inter)}')

    # compare mapping results
    out_log.append(lines='\nCompare mapping results:')

    # both mapped reads
    mapped_wd = set(paf['Woodynano'].keys())
    mapped_py = set(paf['Pychopper'].keys())
    mapped_inter = mapped_wd & mapped_py
    foo = len(mapped_inter)

    # mapped to same loci
    for rname in list(mapped_inter):
        if not paf['Woodynano'][rname].same_loci(paf['Pychopper'][rname]):
            mapped_inter.remove(rname)

    # property of reads mapped to same loci

    out_log.append(f'Both mapped reads: {foo}')
    out_log.append('Mapped to same loci: {same} ({rate})'.format(
        same=len(mapped_inter),
        rate=len(mapped_inter)/foo)
    )
    out_log.append('Woodynano\tPychopper\tStatitics')
    out_log.append(
        '{w}\t{p}\tAverage read length'.format(
            w=np.mean([paf['Woodynano'][i].query_len for i in mapped_inter]),
            p=np.mean([paf['Pychopper'][i].query_len for i in mapped_inter])
            )
        )
    out_log.append(
        '{w}\t{p}\tAverage effective read length'.format(
            w=np.mean([paf['Woodynano'][i].effective_read_length() for i in mapped_inter]),
            p=np.mean([paf['Pychopper'][i].effective_read_length() for i in mapped_inter])
            )
        )
    out_log.append(
        '{w}\t{p}\tAverage residual matched'.format(
            w=np.mean([paf['Woodynano'][i].residual_matched for i in mapped_inter]),
            p=np.mean([paf['Pychopper'][i].residual_matched for i in mapped_inter])
            )
        )
    out_log.append(
        '{w}\t{p}\tAverage align block length'.format(
            w=np.mean([paf['Woodynano'][i].align_block_len for i in mapped_inter]),
            p=np.mean([paf['Pychopper'][i].align_block_len for i in mapped_inter])
            )
        )
    out_log.append(
        '{w}\t{p}\tAverage clipped bases'.format(
            w=np.mean([paf['Woodynano'][i].clipped_base() for i in mapped_inter]),
            p=np.mean([paf['Pychopper'][i].clipped_base() for i in mapped_inter])
            )
        )

    return


def get_stats(species, bio, dirs, data):
    paf = {}
    for software in ['Woodynano', 'Pychopper']:
        # read paf
        # get file path
        foo = get_path(parent=f"{dirs['paf']}{software}/",
                       pattern=f'{species}_{bio}'
                       )
        tmp = PAF.Import(fname=foo).get_primary_paf()
        paf[software] = {get_readname(i): i for i in tmp.alignment}

    # both mapped reads
    mapped_wd = set(paf['Woodynano'].keys())
    mapped_py = set(paf['Pychopper'].keys())
    mapped_inter = mapped_wd & mapped_py
    foo = len(mapped_inter)

    # mapped to same loci
    for rname in list(mapped_inter):
        if not paf['Woodynano'][rname].same_loci(paf['Pychopper'][rname]):
            mapped_inter.remove(rname)

    data['Woodynano'] += [paf['Woodynano'][i].query_len for i in mapped_inter]
    data['Pychopper'] += [paf['Pychopper'][i].query_len for i in mapped_inter]
    data['Species_bio'] += [f'{species}_{bio}']*len(mapped_inter)
    data['Stats'] += ['Read length']*len(mapped_inter)
    data['Readname'] += list(mapped_inter)

    
    data['Woodynano'] += [paf['Woodynano'][i].effective_read_length()
                          for i in mapped_inter]
    data['Pychopper'] += [paf['Pychopper'][i].effective_read_length()
                          for i in mapped_inter]
    data['Species_bio'] += [f'{species}_{bio}']*len(mapped_inter)
    data['Stats'] += ['Effective read length']*len(mapped_inter)
    data['Readname'] += list(mapped_inter)


    data['Woodynano'] += [paf['Woodynano']
                          [i].residual_matched for i in mapped_inter]
    data['Pychopper'] += [paf['Pychopper']
                          [i].residual_matched for i in mapped_inter]
    data['Species_bio'] += [f'{species}_{bio}']*len(mapped_inter)
    data['Stats'] += ['Numbers of matching bases']*len(mapped_inter)
    data['Readname'] += list(mapped_inter)


    data['Woodynano'] += [paf['Woodynano']
                          [i].align_block_len for i in mapped_inter]
    data['Pychopper'] += [paf['Pychopper']
                          [i].align_block_len for i in mapped_inter]
    data['Species_bio'] += [f'{species}_{bio}']*len(mapped_inter)
    data['Stats'] += ['Mapped length (Included gaps)']*len(mapped_inter)
    data['Readname'] += list(mapped_inter)


    data['Woodynano'] += [paf['Woodynano'][i].clipped_base()
                          for i in mapped_inter]
    data['Pychopper'] += [paf['Pychopper'][i].clipped_base()
                          for i in mapped_inter]
    data['Species_bio'] += [f'{species}_{bio}']*len(mapped_inter)
    data['Stats'] += ['Soft clipped bases']*len(mapped_inter)
    data['Readname'] += list(mapped_inter)


    return


def Figure8(data, path, title,lims=[0,6000]):
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
                ylim=lims)

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
                           binwidth=2)
    snsplt2.set(xlabel='WoodyNano-Pychopper', xlim=[-75, 75], ylim=[0,55000])


    fig.savefig(path)

    return


def Figure8_modified(data, path, title, lims=[0, 6000]):
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

#%%
"""
Statistics for paf files
"""
# for species in ['Egr', 'Ptr', 'Lch']:
#     for bio in ['bio1', 'bio2']:
#         compare(species=species, bio=bio, dirs=dirs)

#%%
"""
Plot paf statistics
"""
# data = {'Woodynano': [],
#         'Pychopper': [],
#         'Species_bio': [],
#         'Stats': [],
#         'Readname':[]
#         }
# for species in ['Egr', 'Ptr', 'Lch']:
#     for bio in ['bio1', 'bio2']:
#         get_stats(species=species, bio=bio, dirs=dirs, data=data)

# data = pd.DataFrame(data)
# data.to_csv(
#     '/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figure_stats/Figure8.csv', index=False)

#%%
data = pd.read_csv(
    '/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figure_stats/Figure8.csv'
    )
stats = 'Read length'

data2 = data[data['Stats']==stats]
Figure8_modified(data=data2,
       path=f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figures/Figure8_{stats}.png',
       title=stats)

# # for stats in set(data['Stats']):
# #     Figure8(data=data[data['Stats']==stats], 
# #             path=f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figures/Figure8_{stats}.png',
# #             title=stats)
# stats = 'Mapped length (Included gaps)'

# data2 = data[data['Stats']==stats]
# data2 = data.loc[,Stats] = stats
# data2['Woodynano'] = np.log10(data2['Woodynano'])
# data2['Pychopper'] = np.log10(data2['Pychopper'])

# Figure8(data=data2, 
#        path=f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figures/Figure8_{stats}.png',
#        title=stats,lims=[0,6])



"""
 calculate percentage in the interval
"""
#%%
# from WoodyNano import samtools
# #%%
# importlib.reload(samtools)
# data = pd.read_csv(
#     '/home/woodformation/Processing_data/CCC/WoodyNano_raw/Figure_stats/Figure8.csv'
# )
# data['Diff'] = data['Woodynano']-data['Pychopper']

# data_rl = data[data['Stats'] == 'Read length']
# lower_bound = min(data_rl['Diff'])
# upper_bound = -lower_bound


# between = np.sum((data_rl['Diff'] >= lower_bound) &
#                  (data_rl['Diff'] < upper_bound))/len(data_rl['Diff'])
# between_down = np.sum((data_rl['Diff'] >= lower_bound) &
#                       (data_rl['Diff'] < 0))/len(data_rl['Diff'])
# between_down_bool = np.array((data_rl['Diff'] >= lower_bound) & (data_rl['Diff'] < 0))
# equal = np.sum((data_rl['Diff'] == 0) / len(data_rl['Diff']))
# between_up = np.sum((data_rl['Diff'] > 0) &
#                       (data_rl['Diff'] < upper_bound))/len(data_rl['Diff'])
# over = np.sum((data_rl['Diff'] >= upper_bound) /len(data_rl['Diff']))


# print(f'Between_down: {between_down}\nEqual: {equal}\nBetween_up: {between_up}\nOver: {over}\n')

# """
# Compare soft clip
# """
# data_sf = data[data['Stats']=='Soft clipped bases']
# data_sf_lose = data_sf[np.array(data_sf['Diff']>0)]
# softclip_lose = np.sum(
#     data_sf['Diff'][between_down_bool] > 0) / np.sum(between_down_bool)
# softclip_lose_bool = np.array(data_sf['Diff'] > 0)
# softclip_same = np.sum(data_sf['Diff'][between_down_bool] == 0) / np.sum(between_down_bool)
# softclip_win = np.sum(
#     data_sf['Diff'][between_down_bool] < 0) / np.sum(between_down_bool)
# print(
#     f'softclip_lose: {softclip_lose}\nsoftclip_same: {softclip_same}\nsoftclip_win: {softclip_win}')
# weird_read = data_sf[np.array(
#     [x and y for x, y in zip(softclip_lose_bool, between_down_bool)])]
# weird_read.to_csv('/home/woodformation/Processing_data/CCC/WoodyNano_raw/Weird_read.csv')
# weird_read =  data_sf[np.array([x and y for x, y in zip(softclip_lose_bool, between_down_bool)])][['Species_bio','Readname']]

# for spe in set(weird_read['Species_bio']):
#     tmp = weird_read[weird_read['Species_bio']==spe]['Readname']
    
#     sam_w = samtools.SAM.Import(
#         f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/SAM/Woodynano/{spe}_full_length_woodynano.sam')
    
#     sam_p = samtools.SAM.Import(
#         f'/home/woodformation/Processing_data/CCC/WoodyNano_raw/SAM/Pychopper/{spe}_full_length_pychopper.sam')
    
#     alignment = sam_w.alignment
#     sam_w.alignment = []
#     for a in alignment:
#         if a.qname.split('|')[-1] in list(tmp):
#             sam_w.alignment.append(a)
#     sam_w.export(
#         f'/home/woodformation/Processing_data/CCC/For_EJ/New_SAM/Woodynano/{spe}_woodynano_for_EJ.sam')

#     alignment = sam_p.alignment
#     sam_p.alignment = []
#     for a in alignment:
#         if a.qname.split('|')[-1] in list(tmp):
#             sam_p.alignment.append(a)
#     sam_p.export(
#         f'/home/woodformation/Processing_data/CCC/For_EJ/New_SAM/Pychopper/{spe}_pychopper_for_EJ.sam')


#%%
# out_log = logger(
#     path='/home/woodformation/Processing_data/CCC/WoodyNano_raw/Stats/Difference.tsv')

# out_log.new()
# out_log.append('less\tbetween\tgreater\tstats')
# for stats in set(data['Stats']):
#     tmp = np.array(data[data['Stats'] == stats]['Diff'])
#     out_log.append('{less}\t{between}\t{greater}\t{stats}'.format(less=np.sum(tmp<-50)/len(tmp),
#                                                                   between=np.sum((tmp > -50) & (tmp < 50))/len(tmp),
#                                                                   greater=np.sum(tmp>50)/len(tmp),
#                                                                   stats=stats
#                                                                  )
#                   )



# %%
