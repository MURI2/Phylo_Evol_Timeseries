from __future__ import division
import os
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import phylo_tools as pt
import scipy.stats as stats


def temporal_coverage():
    df = pd.read_csv(pt.get_path() + '/data/bacillus_coverage.txt', sep = '\t', header = 'infer')#, index_col = 0)
    df['cov_ratio'] = df['CP020103'] / df['CP020102']
    df = df.sort_values('Time')
    strains = ['B', 'S']
    treatments = [0,1,2]

    fig = plt.figure()
    count = 0
    #taxa_to_analyze = taxa_to_analyze[:2]
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for treatment in treatments:
        col = pt.get_colors()[str(treatment)]
        for strain in strains:

            ax = fig.add_subplot(3, 2, count+1)
            count+=1

            if (count==1) :
                ax.title.set_text(r'$\mathit{B. \. subtilis} \; \mathrm{sp.} \, \mathrm{168}$')
            if (count==2) :
                ax.title.set_text(r'$\mathit{B. \. subtilis} \; \mathrm{sp.} \, \mathrm{168} \, \Delta \mathrm{spo0A}$')

            reps = list(set(df.Replicate.to_list()))
            for rep in reps:
                df_i = df.loc[(df['Strain'] == strain) & (df['Treatment'] == treatment) & (df['Replicate'] == rep)]
                ax.plot(df_i.Time, df_i.cov_ratio,  'o-', c=col)

            ax.set_ylim([-0.3, 3.5])

            ax.axhline(y=1, color='darkgrey', linestyle='--')

    fig.text(0.5, 0.02, 'Time (days)', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Plasmid-chromosome coverage ratio', va='center', rotation='vertical', fontsize=14)

    fig_name = pt.get_path() + '/figs/plasmid_coverage_plot_new.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_plasmid_coverage():
    df = pd.read_csv(pt.get_path() + '/data/bacillus_coverage.txt', sep = '\t', header = 'infer')#, index_col = 0)
    df['cov_ratio'] = df['CP020103'] / df['CP020102']
    df_group = df.groupby(['Strain', 'Treatment', 'Time'])[["cov_ratio"]].agg(['mean', 'std', 'count'])

    df_0B = df_group.loc[('B', 0), :]
    df_1B = df_group.loc[('B', 1), :]
    df_2B = df_group.loc[('B', 2), :]
    df_0S = df_group.loc[('S', 0), :]
    df_1S = df_group.loc[('S', 1), :]
    df_2S = df_group.loc[('S', 2), :]

    df_0B_time = df_0B.reset_index().Time.values
    df_0B_mean = df_0B.reset_index().cov_ratio['mean'].values
    df_0B_se = df_0B.reset_index().cov_ratio['std'].values / np.sqrt(df_0B.reset_index().cov_ratio['count'].values)

    df_1B_time = df_1B.reset_index().Time.values
    df_1B_mean = df_1B.reset_index().cov_ratio['mean'].values
    df_1B_se = df_1B.reset_index().cov_ratio['std'].values / np.sqrt(df_1B.reset_index().cov_ratio['count'].values)

    df_2B_time = df_2B.reset_index().Time.values
    df_2B_mean = df_2B.reset_index().cov_ratio['mean'].values
    df_2B_se = df_2B.reset_index().cov_ratio['std'].values / np.sqrt(df_2B.reset_index().cov_ratio['count'].values)

    df_0S_time = df_0S.reset_index().Time.values
    df_0S_mean = df_0S.reset_index().cov_ratio['mean'].values
    df_0S_se = df_0S.reset_index().cov_ratio['std'].values / np.sqrt(df_0S.reset_index().cov_ratio['count'].values)

    df_1S_time = df_1S.reset_index().Time.values
    df_1S_mean = df_1S.reset_index().cov_ratio['mean'].values
    df_1S_se = df_1S.reset_index().cov_ratio['std'].values / np.sqrt(df_1S.reset_index().cov_ratio['count'].values)

    df_2S_time = df_2S.reset_index().Time.values
    df_2S_mean = df_2S.reset_index().cov_ratio['mean'].values
    df_2S_se = df_2S.reset_index().cov_ratio['std'].values / np.sqrt(df_2S.reset_index().cov_ratio['count'].values)


    fig = plt.figure()

    plt.subplot(311)
    plt.text(575, 2, "1-day", fontsize=14)
    plt.axhline(y=1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    plt.errorbar(df_0S_time, df_0S_mean, yerr=df_0S_se, fmt='o', label=r'$\Delta spo0A$', mfc='white', color=pt.get_colors()[str(0)], zorder=2)
    plt.errorbar(df_0B_time, df_0B_mean, yerr=df_0B_se, fmt='o', label=r'$wt$', color=pt.get_colors()[str(0)], zorder=3)
    plt.xlim([0,710])
    plt.ylim([0,2.5])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(312)
    plt.text(575, 2, "10-days", fontsize=14)
    plt.axhline(y=1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    plt.errorbar(df_1S_time, df_1S_mean, yerr=df_1S_se, fmt='o', label=r'$\Delta spo0A$', mfc='white', color=pt.get_colors()[str(1)], zorder=2)
    plt.errorbar(df_1B_time, df_1B_mean, yerr=df_1B_se, fmt='o', label=r'$wt$', color=pt.get_colors()[str(1)], zorder=3)
    plt.xlim([0,710])
    plt.ylim([0,2.5])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(313)
    plt.text(575, 2, "100-days", fontsize=14)
    plt.axhline(y=1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    plt.errorbar(df_2S_time, df_2S_mean, yerr=df_2S_se, fmt='o', label=r'$\Delta spo0A$', mfc='white', color=pt.get_colors()[str(2)], zorder=2)
    plt.errorbar(df_2B_time, df_2B_mean, yerr=df_2B_se, fmt='o', label=r'$wt$', color=pt.get_colors()[str(2)], zorder=3)
    plt.xlim([0,710])
    plt.ylim([0,2.5])

    # Set common labels
    fig.text(0.5, 0.04, 'Time (days)', ha='center', va='center', fontsize=18)
    fig.text(0.06, 0.5, 'Plasmid / chromosome coverage', ha='center', va='center', rotation='vertical', fontsize=16)

    fig_name = pt.get_path() + '/figs/plasmid_coverage.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def plot_bPTR():
    df = pd.read_csv(pt.get_path() + '/data/bPTR_clean.txt', sep = '\t', header = 'infer')#, index_col = 0)
    df = df[np.isfinite(df['bPTR'])]

    df_mean = df.groupby(['Strain', 'Treatment', 'Time']).mean().reset_index()

    df_std = df.groupby(['Strain', 'Treatment', 'Time']).agg(stats.variation).reset_index()


    df_0B = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0)]
    df_1B = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 1)]
    df_2B = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 2)]
    df_0S = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0)]
    df_1S = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 1)]
    df_2S = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 2)]


    df_mean_0B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 0)]
    df_mean_1B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 1)]
    df_mean_2B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 2)]
    df_mean_0S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 0)]
    df_mean_1S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 1)]
    df_mean_2S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 2)]


    df_std_0B = df_std.loc[(df_std['Strain'] == 'B') & (df_std['Treatment'] == 0)]
    df_std_1B = df_std.loc[(df_std['Strain'] == 'B') & (df_std['Treatment'] == 1)]
    df_std_2B = df_std.loc[(df_std['Strain'] == 'B') & (df_std['Treatment'] == 2)]
    df_std_0S = df_std.loc[(df_std['Strain'] == 'S') & (df_std['Treatment'] == 0)]
    df_std_1S = df_std.loc[(df_std['Strain'] == 'S') & (df_std['Treatment'] == 1)]
    df_std_2S = df_std.loc[(df_std['Strain'] == 'S') & (df_std['Treatment'] == 2)]


    fig = plt.figure()

    plt.subplot(311)
    plt.axhline( y=np.mean(df_std_0B.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_std_0S.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_std_0B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_std_0S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 0.65, "1-day", fontsize=14)
    plt.scatter(df_std_0B.Time.values, df_std_0B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(0)], zorder=3)
    plt.scatter(df_std_0S.Time.values, df_std_0S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(0)], zorder=4)
    plt.xlim([0,710])
    plt.ylim([0,0.8])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(312)
    plt.axhline( y=np.mean(df_std_1B.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_std_1S.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_std_1B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_std_1S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 0.65, "10-days", fontsize=14)
    plt.scatter(df_std_1B.Time.values, df_std_1B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(1)], zorder=3)
    plt.scatter(df_std_1S.Time.values, df_std_1S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(1)], zorder=4)
    plt.xlim([0,710])
    plt.ylim([0,0.8])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(313)
    plt.axhline( y=np.mean(df_std_2B.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_std_2S.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_std_2B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_std_2S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 0.65, "100-days", fontsize=14)
    plt.scatter(df_std_2B.Time.values, df_std_2B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(2)])
    plt.scatter(df_std_2S.Time.values, df_std_2S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(2)])
    plt.xlim([0,710])
    plt.ylim([0,0.8])

    # Set common labels
    fig.text(0.5, 0.04, 'Time (days)', ha='center', va='center', fontsize=18)
    fig.text(0.06, 0.5, r'$\frac{\sigma}{\mu}$' + ' peak-to-trough ratio', ha='center', va='center', rotation='vertical', fontsize=18)

    fig_name = pt.get_path() + '/figs/bptr_cv_bacillus.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




    fig = plt.figure()

    plt.subplot(311)
    plt.axhline( y=np.mean(df_mean_0B.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_0S.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_0B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_0S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1.65, "1-day", fontsize=14)
    plt.scatter(df_mean_0B.Time.values, df_mean_0B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(0)], zorder=3)
    plt.scatter(df_mean_0S.Time.values, df_mean_0S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(0)], zorder=4)
    plt.xlim([0,710])
    plt.ylim([0.2,2])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(312)
    plt.axhline( y=np.mean(df_mean_1B.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_1S.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_1B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_1S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1.65, "10-days", fontsize=14)
    plt.scatter(df_mean_1B.Time.values, df_mean_1B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(1)], zorder=3)
    plt.scatter(df_mean_1S.Time.values, df_mean_1S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(1)], zorder=4)
    plt.xlim([0,710])
    plt.ylim([0.4,2.5])
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.subplot(313)
    plt.axhline( y=np.mean(df_mean_2B.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_2S.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_2B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_2S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1, "100-days", fontsize=14)
    plt.scatter(df_mean_2B.Time.values, df_mean_2B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(2)])
    plt.scatter(df_mean_2S.Time.values, df_mean_2S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(2)])
    plt.xlim([0,710])
    plt.ylim([0.5,1.2])

    # Set common labels
    fig.text(0.5, 0.04, 'Time (days)', ha='center', va='center', fontsize=18)
    fig.text(0.06, 0.5, 'Mean peak-to-trough ratio', ha='center', va='center', rotation='vertical', fontsize=18)

    fig_name = pt.get_path() + '/figs/bptr_mean_bacillus.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_bPTR_all():
    df = pd.read_csv(pt.get_path() + '/data/bPTR_clean.txt', sep = '\t', header = 'infer')#, index_col = 0)
    df = df[np.isfinite(df['bPTR'])]

    df_mean = df.groupby(['Strain', 'Treatment', 'Time']).mean().reset_index()

    df_std = df.groupby(['Strain', 'Treatment', 'Time']).agg(stats.variation).reset_index()

    df_mean_0B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 0)]
    df_mean_1B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 1)]
    df_mean_2B = df_mean.loc[(df_mean['Strain'] == 'B') & (df_mean['Treatment'] == 2)]
    df_mean_0S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 0)]
    df_mean_1S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 1)]
    df_mean_2S = df_mean.loc[(df_mean['Strain'] == 'S') & (df_mean['Treatment'] == 2)]


    fig = plt.figure()

    plt.axhline( y=np.mean(df_mean_0B.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_0S.bPTR.values), color=pt.get_colors()[str(0)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_0B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_0S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1.65, "1-day", fontsize=14)
    plt.scatter(df_mean_0B.Time.values, df_mean_0B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(0)], zorder=3)
    plt.scatter(df_mean_0S.Time.values, df_mean_0S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(0)], zorder=4)

    plt.axhline( y=np.mean(df_mean_1B.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_1S.bPTR.values), color=pt.get_colors()[str(1)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_1B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_1S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1.65, "10-days", fontsize=14)
    plt.scatter(df_mean_1B.Time.values, df_mean_1B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(1)], zorder=3)
    plt.scatter(df_mean_1S.Time.values, df_mean_1S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(1)], zorder=4)

    plt.axhline( y=np.mean(df_mean_2B.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=1)
    plt.axhline( y=np.mean(df_mean_2S.bPTR.values), color=pt.get_colors()[str(2)], linestyle=':', alpha = 0.8, zorder=2)
    plt.text(5, np.mean(df_mean_2B.bPTR.values), r'$\mathrm{wt}$')
    plt.text(5, np.mean(df_mean_2S.bPTR.values) + 0.05, r'$\Delta\mathrm{spo0A}$')
    plt.text(575, 1, "100-days", fontsize=14)
    plt.scatter(df_mean_2B.Time.values, df_mean_2B.bPTR.values, label=r'$wt$', color=pt.get_colors()[str(2)])
    plt.scatter(df_mean_2S.Time.values, df_mean_2S.bPTR.values, label=r'$\Delta spo0A$', facecolors='none', color=pt.get_colors()[str(2)])

    #plt.xlim([0,710])
    #plt.ylim([0.5,1.2])

    # Set common labels
    plt.xlabel( 'Time (days)', fontsize=18)
    plt.ylabel('Mean peak-to-trough ratio', fontsize=18)

    fig_name = pt.get_path() + '/figs/bptr_mean_bacillus_all.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_allele_freqs():



#plot_bPTR_all()
#temporal_coverage()
