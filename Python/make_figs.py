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


    df_0B1 = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0) & (df['Replicate'] == 1)]
    df_0B2 = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0) & (df['Replicate'] == 2)]
    df_0B3 = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0) & (df['Replicate'] == 3)]
    df_0B4 = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0) & (df['Replicate'] == 4)]
    df_0B5 = df.loc[(df['Strain'] == 'B') & (df['Treatment'] == 0) & (df['Replicate'] == 5)]

    df_0S1 = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0) & (df['Replicate'] == 1)]
    df_0S2 = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0) & (df['Replicate'] == 2)]
    df_0S3 = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0) & (df['Replicate'] == 3)]
    df_0S4 = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0) & (df['Replicate'] == 4)]
    df_0S5 = df.loc[(df['Strain'] == 'S') & (df['Treatment'] == 0) & (df['Replicate'] == 5)]

    fig = plt.figure()

    plt.plot(df_0B1.Time, df_0B1.cov_ratio, 'o-')
    plt.plot(df_0B2.Time, df_0B2.cov_ratio, 'o-')
    plt.plot(df_0B3.Time, df_0B3.cov_ratio, 'o-')
    plt.plot(df_0B4.Time, df_0B4.cov_ratio, 'o-')
    plt.plot(df_0B5.Time, df_0B5.cov_ratio, 'o-')

    plt.xlabel('Time (days)', fontsize = 18)
    plt.ylabel('Plasmid / chromosome coverage', fontsize = 16)

    fig_name = pt.get_path() + '/figs/plasmid_coverage_plot.png'
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






#plot_bPTR()
temporal_coverage()
