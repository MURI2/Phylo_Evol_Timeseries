from __future__ import division
import os, sys, json
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import phylo_tools as pt
import scipy.stats as stats
import statsmodels.api as sm

import parse_file
import timecourse_utils
import mutation_spectrum_utils

from scipy.special import gammaln

from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

np.random.seed(123456789)

treatments=pt.treatments
replicates = pt.replicates


latex_formal_dict = {  'B': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610}$',
                'S': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610} \, \Delta \mathrm{spo0A} $',
                'C': r'$\mathit{Caulobacter \, crescentus} \; \mathrm{NA1000}$',
                'D': r'$\mathit{Deinococcus \, radiodurans} \; \mathrm{BAA-816}$',
                'P': r'$\mathit{Pseudomonas \,} \; \mathrm{sp. \, KBS0710}$',
                'F': r'$\mathit{Pedobacter \,} \; \mathrm{sp. \, KBS0701}$',
                'J': r'$\mathit{Janthinobacterium \,} \; \mathrm{sp. \, KBS0711}$'
                }



latex_dict = {  'B': r'$\mathit{Bacillus\, subtilis} \, \mathrm{wt} $',
                'S': r'$\mathit{Bacillus\, subtilis} \, \Delta \mathrm{spo0A} $',
                'C': r'$\mathit{Caulobacter \, crescentus}$',
                'D': r'$\mathit{Deinococcus \, radiodurans}$',
                'P': r'$\mathit{Pseudomonas \,} \; \mathrm{sp. \, KBS0710}$',
                'F': r'$\mathit{Pedobacter \,} \; \mathrm{sp. \, KBS0701}$',
                'J': r'$\mathit{Janthinobacterium \,} \; \mathrm{sp. \, KBS0711}$'
                }





def plot_spores():
    path_IN = pt.get_path() + '/data/spore_assay/Sporulation_170912_long.txt'
    IN = pd.read_csv(path_IN, sep = '\t')
    IN = IN.loc[IN['Time_hours'] <= 400]
    #d100
    IN_0B1_100 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 100)]
    IN_2B1_100 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 100)]
    IN_mean_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).mean().reset_index()
    IN_mean_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).mean().reset_index()
    IN_std_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).std().reset_index()
    IN_std_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).std().reset_index()
    # Day 500
    IN_0B1_500 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 500)]
    IN_2B1_500 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 500)]
    IN_mean_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).mean().reset_index()
    IN_mean_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).mean().reset_index()
    IN_std_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).std().reset_index()
    IN_std_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).std().reset_index()

    fig = plt.figure(figsize=(6,3))

    #plt.scatter(IN_mean_0B1.Time_hours.values, IN_mean_0B1.Vegetative_percent.values, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
    plt.plot(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
        'b-',  c = '#FF6347')
    plt.errorbar(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
        IN_std_0B1_100.Vegetative_percent.values,  linestyle='None', marker='o', c='#87CEEB', elinewidth=1.5, label="1-day WT, day 100",)
    plt.errorbar(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
        IN_std_2B1_100.Vegetative_percent.values, linestyle='None', marker='o', c = '#FF6347', elinewidth=1.5, label="100-day WT, day 100",)

    plt.plot(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        'b-',  c = '#FF6347')
    plt.errorbar(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        IN_std_0B1_500.Vegetative_percent.values,  linestyle='None', marker='v', c='#87CEEB', elinewidth=1.5, label="1-day WT, day 500",)
    plt.errorbar(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        IN_std_2B1_500.Vegetative_percent.values, linestyle='None', marker='v', c = '#FF6347', elinewidth=1.5, label="100-day WT, day 500",)
    #plt.title('Bacillus sporulation', fontsize = 24)
    plt.xlabel('Time (hours)', fontsize = 18)
    plt.ylabel('Percent spores', fontsize = 18)
    plt.ylim(0.0008, 1.1)
    plt.xlim(-5, 400)
    plt.yscale('log')
    plt.legend(numpoints=1, prop={'size':8},  loc='lower right', frameon=False)

    plt.title(latex_dict['B'], fontsize=24)

    fig_name = pt.get_path() + '/figs/spore_assay.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_spo0a_fitness():
    df = pd.read_csv(pt.get_path() + '/data/competition_2018-1-9-count.txt', sep = '\t')
    rows_to_keep = []
    for index, row in df.iterrows():
        wt = row['WT']
        spo0a = row['spoA']
        if (wt == 'TMTC') or (spo0a == 'TMTC') :
            continue
        wt_spo0a = int(wt) + int(spo0a)
        if (wt_spo0a < 30):
            continue
        rows_to_keep.append(index)


    def cfus_ml(column, conc):
        return column * (10 ** (int(conc) * -1))

    df_keep = df.ix[rows_to_keep]
    df_keep.WT = df_keep.WT.astype(int)
    df_keep.spoA = df_keep.spoA.astype(int)

    df_keep['WT_cfus_ml'] = df_keep.apply(lambda x: cfus_ml(x.WT, x.Concentration), axis=1)
    df_keep['spoA_cfus_ml'] = df_keep.apply(lambda x: cfus_ml(x.spoA, x.Concentration), axis=1)

    df_keep = df_keep.drop(['Concentration', 'WT', 'spoA', 'Rep'], 1)
    df_keep = df_keep.groupby(['Day','Flask'], as_index=False).mean()

    flask_1 = df_keep.loc[df_keep['Flask'] == 1]
    flask_2 = df_keep.loc[df_keep['Flask'] == 2]
    flask_3 = df_keep.loc[df_keep['Flask'] == 3]

    relative_fitness_1 = np.log((flask_1['spoA_cfus_ml'].values / flask_1['WT_cfus_ml'].values) * ( 0.51/0.49))
    relative_fitness_2 = np.log((flask_2['spoA_cfus_ml'].values / flask_2['WT_cfus_ml'].values) * ( 0.48/0.52))
    relative_fitness_3 = np.log((flask_3['spoA_cfus_ml'].values / flask_3['WT_cfus_ml'].values) * ( 0.54/0.46))

    relative_fitness_per_time_1 = np.log((flask_1['spoA_cfus_ml'].values / flask_1['WT_cfus_ml'].values) * ( 0.51/0.49)) /  flask_1['Day'].values
    relative_fitness_per_time_2 = np.log((flask_2['spoA_cfus_ml'].values / flask_2['WT_cfus_ml'].values) * ( 0.48/0.52)) /  flask_2['Day'].values
    relative_fitness_per_time_3 = np.log((flask_3['spoA_cfus_ml'].values / flask_3['WT_cfus_ml'].values) * ( 0.54/0.46)) /  flask_3['Day'].values


    zipped_relative = list(zip(list(relative_fitness_1), list(relative_fitness_2), list(relative_fitness_3)))
    relative_mean = (relative_fitness_1 + relative_fitness_2 + relative_fitness_3) / 3
    relative_se_list = []
    for i , item in enumerate(zipped_relative):
        relative_se_list.append(2*np.std(np.asarray(item)) / np.sqrt(len(item)))

    zipped_relative_time = list(zip(list(relative_fitness_per_time_1), list(relative_fitness_per_time_2), list(relative_fitness_per_time_3)))
    relative_time_mean = (relative_fitness_per_time_1 + relative_fitness_per_time_2 + relative_fitness_per_time_3) / 3
    relative_time_se_list = []
    for i , item in enumerate(zipped_relative_time):
        relative_time_se_list.append(2*np.std(np.asarray(item)) / np.sqrt(len(item)))


    fig = plt.figure(figsize = (10, 9))

    ax_relative = plt.subplot2grid((2, 1), (0, 0), colspan=1)
    ax_time_relative = plt.subplot2grid((2, 1), (1, 0), colspan=1)

    ax_relative.axhline(y=0, color='grey', linestyle='--', lw = 3)
    ax_relative.errorbar(flask_1['Day'].values, relative_mean, relative_se_list, linestyle='-', marker='o', lw = 3)
    ax_relative.set_ylim(-5, 5)
    ax_relative.set_ylabel(latex_dict['S'] + '\nrelative fitness at ' + r'$t$' + ', ' + r'$X(t)$'  , fontsize = 16)
    ax_relative.set_xscale('log', basex=10)

    ax_time_relative.axhline(y=0, color='grey', linestyle='--', lw = 3)
    ax_time_relative.errorbar(flask_1['Day'].values, relative_time_mean, relative_time_se_list, linestyle='-', marker='o', lw = 3)
    ax_time_relative.set_ylim(-0.35, 0.75)
    ax_time_relative.set_ylabel(latex_dict['S'] + '\nrelative fitness, ' + r'$\Delta X$'  , fontsize = 16)
    ax_time_relative.set_xscale('log', basex=10)


    ax_time_relative.set_xlabel('Days, ' + r'$t$', fontsize = 20)

    fig_name = pt.get_path() + '/figs/fitness_spo0a.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def temporal_coverage_B_S():

    taxa = ['B', 'S']

    fig = plt.figure(figsize = (14, 12))

    for taxon_idx, taxon in enumerate(taxa):

        for treatment_idx, treatment in enumerate(treatments):

            ax_i = plt.subplot2grid((3, 2), (treatment_idx, taxon_idx), colspan=1)

            if treatment_idx == 0:
                ax_i.set_title(latex_dict[taxon] + '\nplasmid pBS32 (84,215bp)', fontsize=20 )

            if taxon_idx == 0:
                ax_i.set_ylabel( str(int(10** int(treatment) ) ) + '-day transfer', fontsize=20  )


            for replicate in replicates:

                population = '%s%s%s' % (treatment, taxon, replicate)

                coverage_ratio = []
                times = []

                for time in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:

                    sample = '%s_%s' % (population, time)
                    sample_path = pt.get_path() + '/data/rebreseq_json/%s.json' % sample

                    if os.path.isfile(sample_path) == False:
                        continue

                    data = json.load(open(sample_path))

                    if data['references']['reference']['NZ_CP020102']['coverage_average'] >= 50:
                        coverage_ratio.append(data['references']['reference']['NZ_CP020103']['coverage_average'] / data['references']['reference']['NZ_CP020102']['coverage_average'])
                        times.append(time)

                ax_i.plot(times, coverage_ratio,  marker=pt.plot_species_marker(taxon), c=pt.get_colors(treatment), linestyle='dashed', alpha=0.8, ms = 10, lw=1.5, zorder=2)
                ax_i.axhline(y=1, color='darkgrey', linestyle='--', zorder=1)
                ax_i.axhline(y=0, color='k', linestyle=':', zorder=1)
                ax_i.set_ylim([-0.3, 3.5])
                if taxon == 'B':
                    ax_i.xaxis.set_tick_params(labelsize=10)

    fig.text(0.5, 0.05, 'Days', ha='center', fontsize=30)
    fig.text(0.05, 0.4, 'Plasmid-chromosome coverage ratio', va='center', rotation='vertical', fontsize=28)

    fig_name = pt.get_path() + '/figs/plasmid_coverage_B_S.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






def temporal_coverage_D():

    plasmids = ['NC_001264', 'NC_000958', 'NC_000959']

    pladmid_name_dict = {'NC_001264':'chromosome 2 (412,348bp)', 'NC_000958': 'plasmid MP1 (177,466bp)', 'NC_000959': 'plasmid CP1 (45,704bp)'}

    fig = plt.figure(figsize = (18, 14))

    taxon = 'D'

    for plasmid_idx, plasmid in enumerate(plasmids):

        for treatment_idx, treatment in enumerate(treatments):

            ax_i = plt.subplot2grid((3, 3), (treatment_idx, plasmid_idx), colspan=1)

            if treatment_idx == 0:

                ax_i.set_title( latex_dict[taxon] + '\n' + pladmid_name_dict[plasmid], fontsize=20 )

            if plasmid_idx == 0:

                ax_i.set_ylabel( str(int(10** int(treatment) ) ) + '-day transfer', fontsize=20  )

            for replicate in replicates:

                population = '%s%s%s' % (treatment, taxon, replicate)

                coverage_ratio = []
                times = []

                for time in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:

                    sample = '%s_%s' % (population, str(time))
                    sample_path = pt.get_path() + '/data/rebreseq_json/%s.json' % sample

                    if os.path.isfile(sample_path) == False:
                        continue

                    data = json.load(open(sample_path))

                    if data['references']['reference']['NC_001263']['coverage_average'] >= 50:
                        coverage_ratio.append(data['references']['reference'][plasmid]['coverage_average'] / data['references']['reference']['NC_001263']['coverage_average'])
                        times.append(time)

                ax_i.plot(times, coverage_ratio,  marker=pt.plot_species_marker(taxon), c=pt.get_colors(treatment), linestyle='dashed', alpha=0.8, ms = 10, lw=1.5, zorder=2)
                ax_i.axhline(y=1, color='darkgrey', linestyle='--', zorder=1)
                ax_i.axhline(y=0, color='k', linestyle=':', zorder=1)

                if treatment_idx == plasmid_idx == 2:
                    ax_i.set_ylim([-0.3, 8.5])
                else:
                    ax_i.set_ylim([-0.3, 5.5])


    fig.text(0.5, 0.05, 'Days', ha='center', fontsize=30)
    fig.text(0.05, 0.5, 'Plasmid-chromosome coverage ratio', va='center', rotation='vertical', fontsize=28)

    fig_name = pt.get_path() + '/figs/plasmid_coverage_D.png'
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



def plot_allele_freqs_all_treats(strain):

    fig = plt.figure(figsize = (12, 6))

    row_count = 0
    for treatment in treatments:

        column_count = 0

        for replicate in replicates:
            ax_i = plt.subplot2grid((3, 5), (row_count, column_count), colspan=1)

            population = treatment + strain + replicate
            print(population)
            annotated_timecourse_path = pt.get_path() + "/data/timecourse_final/%s_annotated_timecourse.txt" % population
            if os.path.exists(annotated_timecourse_path) == False:
                continue
            annotated_timecourse_file = open(annotated_timecourse_path ,"r")

            first_line = annotated_timecourse_file.readline()
            first_line = first_line.strip()
            first_line_items = first_line.split(",")
            times = np.asarray([float(x.strip().split(':')[1]) for x in first_line_items[13::2]])
            times = np.insert(times, 0, 0, axis=0)


            for i, line in enumerate(annotated_timecourse_file):
                line = line.strip()
                items = line.split(",")
                pass_or_fail = items[12].strip()
                if pass_or_fail == 'FAIL':
                    continue

                alt_cov = np.asarray([ float(x) for x in items[13::2]])
                total_cov = np.asarray([ float(x) for x in items[14::2]])
                # pseudocount to avoid divide by zero error
                alt_cov = alt_cov
                total_cov = total_cov + 1
                freqs = alt_cov / total_cov
                freqs = np.insert(freqs, 0, 0, axis=0)

                rgb = pt.mut_freq_colormap()
                rgb = pt.lighten_color(rgb, amount=0.5)


                if len(times) == len(freqs) + 1:
                    freqs = np.insert(freqs, 0, 0, axis=0)


                ax_i.plot(times, freqs, '.-', c=rgb, alpha=0.4)


            ax_i.set_xlim([0,max(times)])
            ax_i.set_ylim([0, 1])

            if column_count == 0:
                ax_i.set_ylabel( str(10**int(treatment)) + '-day transfers', fontsize =12  )

            ax_i.tick_params(axis="x", labelsize=8)
            ax_i.tick_params(axis="y", labelsize=8)

            column_count += 1
        row_count +=1


    fig.text(0.5, 0.04, 'Days, ' + r'$t$', ha='center', va='center', fontsize=24)
    fig.text(0.05, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=24)


    fig.suptitle(latex_dict[strain], fontsize=28, fontweight='bold')
    fig_name = pt.get_path() + '/figs/mut_trajectories_%s.png'
    fig.savefig(fig_name % strain, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()








# calculate_unnormalized_survival_from_vector function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_unnormalized_survival_from_vector(xs, min_x=None, max_x=None, min_p=1e-10):
    if min_x==None:
        min_x = xs.min()-1

    if max_x==None:
        max_x = xs.max()+1

    unique_xs = set(xs)
    unique_xs.add(min_x)
    unique_xs.add(max_x)

    xvalues = []
    num_observations = []

    for x in sorted(unique_xs):
        xvalues.append(x)
        num_observations.append( (xs>=x).sum() )

    # So that we can plot CDF, SF on log scale
    num_observations[0] -= min_p
    num_observations[1] -= min_p
    num_observations[-1] += min_p

    return np.array(xvalues), np.array(num_observations)







def get_mutation_fixation_trajectories(population):

    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    times = mutations[0][9]
    Ms = np.zeros_like(times)*1.0
    fixed_Ms = np.zeros_like(times)*1.0

    #transit_times[population] = []

    for mutation_idx in range(0,len(mutations)):

        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

        state_Ls = state_trajectories[mutation_idx]

        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]

        t0,tf,transit_time = timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
        if t0==tf==transit_time==None:
            continue

        #print(masked_times, masked_freqs)

        interpolating_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)

        fs = interpolating_function(times)
        fs[fs<0]=0

        # Record
        Ms += fs
        if masked_state_Ls[-1] in parse_file.well_mixed_fixed_states:
            fixed_Ms += (times>=tf)


    return times, Ms, fixed_Ms





def plot_B_S_mutation_trajectory():
    sys.stderr.write("Loading mutation data...\n")

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    delta_mutation_trajectories = {}
    #transit_times = {}
    taxa = ['B', 'S']

    for treatment in treatments:
        for taxon in taxa:
            for replicate in replicates:

                population = treatment + taxon + replicate
                sys.stderr.write("Processing %s...\t" % population)

                times, Ms, fixed_Ms = get_mutation_fixation_trajectories(population)

                fixed_mutation_trajectories[population] = (times, fixed_Ms)
                mutation_trajectories[population] = (times,np.log10(Ms))
                delta_mutation_trajectories[population] = (times[1:], np.log10(Ms[1:]/Ms[:-1] ))

                sys.stderr.write("analyzed %d mutations!\n" % len(Ms))


    fig = plt.figure(figsize = (10, 9))

    column_count = 0

    for treatment in treatments:

        ax_t_vs_M = plt.subplot2grid((3, 3), (0, column_count), colspan=1)

        ax_t_vs_delta_M = plt.subplot2grid((3, 3), (1, column_count), colspan=1)

        ax_M_vs_F = plt.subplot2grid((3, 3), (2, column_count), colspan=1)

        #ax_M_vs_F.plot([0,300],[0,300],'--',linewidth=1,color='k', zorder=1)

        for taxon_i, taxon in enumerate(taxa):

            treatment_taxon_populations = []

            for replicate in replicates:

                population = treatment + taxon + replicate

                Mts,Ms = mutation_trajectories[population]
                fixed_Mts, fixed_Ms = fixed_mutation_trajectories[population]
                deltaMts, deltaMs = delta_mutation_trajectories[population]

                ax_t_vs_M.plot(Mts, 10**Ms, 'o-',color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
                ax_t_vs_M.set_yscale('log', basey=10)
                ax_t_vs_M.tick_params(axis='x', labelsize=8)

                # back transform to format plot axes
                ax_t_vs_delta_M.plot(deltaMts, 10**deltaMs, color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon))
                ax_t_vs_delta_M.set_yscale('log', basey=10)

                ax_M_vs_F.plot(fixed_Mts, fixed_Ms, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
                #ax_M_vs_F.set_xlabel('Days, ' + r'$t$', fontsize = 12)

                treatment_taxon_populations.append(population)

            avg_Mts, avg_Ms = timecourse_utils.average_trajectories([mutation_trajectories[population] for population in treatment_taxon_populations])

            avg_deltaMts, avg_deltaMs = timecourse_utils.average_trajectories([delta_mutation_trajectories[population] for population in treatment_taxon_populations])


            if taxon == 'B':
                ls = '--'
            else:
                ls = ':'

            ax_t_vs_delta_M.axhline(y=1, c='grey', linestyle=':', lw=3, zorder=1)
            ax_t_vs_M.plot(avg_Mts, 10**avg_Ms, ls,color='k', marker=" ", alpha=1, linewidth=4, zorder=2)
            ax_t_vs_delta_M.plot(avg_deltaMts, 10**avg_deltaMs, ls,color='k', marker=" ", alpha=1, linewidth=4, zorder=2)


            if (taxon_i == 0) and (column_count==0):
                legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}_{WT} (t)$'),
                                   Line2D([0], [0], ls=':', color='k', lw=1.5, label= r'$\overline{M}_{\Delta \mathrm{spo0A}} (t)$')]
                ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=8)

        ax_t_vs_M.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

        if treatment == '2':
            ax_M_vs_F.yaxis.set_major_locator(MaxNLocator(integer=True))

        if column_count == 0:

            ax_t_vs_M.set_ylabel('Mutations, ' + r'$M(t)$', fontsize = 15)
            ax_M_vs_F.set_ylabel('Fixed mutations', fontsize = 15)
            ax_t_vs_delta_M.set_ylabel('Change in mutations,\n' + r'$M(t)/M(t-1)$', fontsize = 15)

        column_count += 1

    fig.text(0.53, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=28)


    fig_name = pt.get_path() + '/figs/rate_B_S.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_taxon_mutation_trajectory(taxon):

    sys.stderr.write("Loading mutation data...\n")

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    delta_mutation_trajectories = {}
    #transit_times = {}

    for treatment in treatments:
        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            sys.stderr.write("Processing %s...\t" % population)

            times, Ms, fixed_Ms = get_mutation_fixation_trajectories(population)

            fixed_mutation_trajectories[population] = (times, fixed_Ms)
            mutation_trajectories[population] = (times,np.log10(Ms))
            delta_mutation_trajectories[population] = (times[1:], np.log10(Ms[1:]/Ms[:-1] ))

            sys.stderr.write("analyzed %d mutations!\n" % len(Ms))

    fig = plt.figure(figsize = (10, 9))

    column_count = 0

    for treatment in treatments:

        ax_t_vs_M = plt.subplot2grid((3, 3), (0, column_count), colspan=1)

        ax_t_vs_delta_M = plt.subplot2grid((3, 3), (1, column_count), colspan=1)

        ax_M_vs_F = plt.subplot2grid((3, 3), (2, column_count), colspan=1)

        treatment_taxon_populations = []

        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            Mts,Ms = mutation_trajectories[population]
            fixed_Mts, fixed_Ms = fixed_mutation_trajectories[population]
            deltaMts, deltaMs = delta_mutation_trajectories[population]

            ax_t_vs_M.plot(Mts, 10**Ms, 'o-',color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            ax_t_vs_M.set_yscale('log', basey=10)
            ax_t_vs_M.tick_params(axis='x', labelsize=8)

            # back transform to format plot axes
            ax_t_vs_delta_M.plot(deltaMts, 10**deltaMs, color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon))
            ax_t_vs_delta_M.set_yscale('log', basey=10)

            ax_M_vs_F.plot(fixed_Mts, fixed_Ms, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            #ax_M_vs_F.set_xlabel('Days, ' + r'$t$', fontsize = 12)

            treatment_taxon_populations.append(population)

        avg_Mts, avg_Ms = timecourse_utils.average_trajectories([mutation_trajectories[population] for population in treatment_taxon_populations])

        avg_deltaMts, avg_deltaMs = timecourse_utils.average_trajectories([delta_mutation_trajectories[population] for population in treatment_taxon_populations])

        ax_t_vs_delta_M.axhline(y=1, c='grey', linestyle=':', lw=3, zorder=1)
        ax_t_vs_M.plot(avg_Mts, 10**avg_Ms, '--',color='k', marker=" ", alpha=1, linewidth=4, zorder=2)
        ax_t_vs_delta_M.plot(avg_deltaMts, 10**avg_deltaMs, '--',color='k', marker=" ", alpha=1, linewidth=4, zorder=2)


        if (column_count==0):
            legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}(t)$')]
            ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=8)

        ax_t_vs_M.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

        #if treatment == '2':
        #    ax_M_vs_F.yaxis.set_major_locator(MaxNLocator(integer=True))

        if column_count == 0:

            ax_t_vs_M.set_ylabel('Mutations, ' + r'$M(t)$', fontsize = 15)
            ax_M_vs_F.set_ylabel('Fixed mutations', fontsize = 15)
            ax_t_vs_delta_M.set_ylabel('Change in mutations,\n' + r'$M(t)/M(t-1)$', fontsize = 15)

        column_count += 1

    fig.text(0.53, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=28)


    fig.suptitle(latex_dict[taxon], fontsize=30)


    fig_name = pt.get_path() + '/figs/rate_%s.png' % taxon
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def allele_survival():
    frequencies = np.linspace(0.1,0.9,201)
    df = 0.05
    fstar = 0.5

    #tstar = 20025

    null_num_in_bin = np.zeros_like(frequencies)
    null_avg_f = np.zeros_like(frequencies)

    origin_fixation_times = {}

    taxa = ['B', 'S']

    #fig = plt.figure()
    fig = plt.figure(figsize = (12, 6))

    for treatment in ['0']:
        for taxon_idx, taxon in enumerate(taxa):
            ax_i = plt.subplot2grid((1, 2), (0, taxon_idx), colspan=1)
            for replicate in ['1', '2', '3', '4', '5']:
                population = treatment + taxon + replicate

                if population in pt.populations_to_ignore:
                    continue

                num_runs = []

                num_in_bin = np.zeros_like(null_num_in_bin)
                avg_f = np.zeros_like(null_avg_f)

                origin_fixation_times[population] = ([],[],[])

                sys.stderr.write("Processing fixation probabilities for %s...\n" % population)

                # load mutations
                mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
                population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
                state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

                num_runs = []

                for mutation_idx in range(0,len(mutations)):

                    location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

                    state_Ls = state_trajectories[mutation_idx]

                    good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                    freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                    masked_times = times[good_idxs]
                    masked_freqs = freqs[good_idxs]
                    masked_depths = depths[good_idxs]
                    masked_state_Ls = state_Ls[good_idxs]

                    # Estimate appearance and fixation times
                    if masked_state_Ls[-1] in parse_file.well_mixed_fixed_states:
                        t0,tf,dt = timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
                        origin_fixation_times[population][0].append(t0)
                        origin_fixation_times[population][1].append(tf)
                        origin_fixation_times[population][2].append(dt)

                    # Now split the trajectory into independent polymorphic runs
                    # each of which contains a single final state (unless end of experiment)
                    independent_runs = timecourse_utils.split_well_mixed_hmm(masked_times,masked_freqs, masked_state_Ls)
                    num_runs.append(len(independent_runs))

                    for run_idxs in independent_runs:

                        if len(run_idxs)<2:
                            # need at least one polymorphic state and one final state
                            continue
                        # initial time
                        t = masked_times[run_idxs[0]]

                        # get final state
                        final_state = masked_state_Ls[run_idxs[-1]]

                        # get frequency of parent clade during run
                        #if final_state == parse_file.well_mixed_hmm_states['F'] or final_state==parse_file.well_mixed_hmm_states['P']:
                        #    parent_freqs = masked_freqs
                        #else:
                        #parent_freqs = masked_freqs
                        #elif final_state == parse_file.clade_hmm_states['F'] or final_state == parse_file.well_mixed_hmm_states['P']:
                        #    parent_freqs = masked_fmajors
                        #elif final_state == parse_file.clade_hmm_states['F'] or final_state == parse_file.clade_hmm_states['Pm']:
                        #    parent_freqs = masked_fmajors
                        #else:
                        #    parent_freqs = masked_fextincts


                        # don't neet to renormalize the freqs because we have no population structure

                        #renormalized_freqs = np.clip(masked_freqs[run_idxs]/parent_freqs[run_idxs],0,1)

                        # get fixation weight
                        if final_state in parse_file.well_mixed_fixed_states:
                            fixation_weight = 1.0
                        elif final_state in parse_file.well_mixed_extinct_states:
                            fixation_weight = 0.0
                        else:
                            fixation_weight = masked_freqs[-1] > fstar

                        individual_bin_weights = np.exp(-np.power((masked_freqs[:,None]-frequencies[None,:])/df,2))
                        individual_bin_weights *= (individual_bin_weights>0.01)

                        bin_weights = individual_bin_weights.sum(axis=0)
                        fixation_weights = (individual_bin_weights*fixation_weight).sum(axis=0)

                        #num_in_bin += bin_weights
                        num_in_bin += bin_weights
                        avg_f += fixation_weights

                        null_num_in_bin += bin_weights
                        null_avg_f += fixation_weights


                avg_f = avg_f/(num_in_bin+(num_in_bin<1))

                ax_i.plot(frequencies[(num_in_bin>=1)], avg_f[(num_in_bin>=1)], '.-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=0.9, linewidth=2.5, markersize = 8, zorder=2)
                ax_i.plot([0, 1], [0, 1], '--', c = 'black', alpha=0.9, markersize = 10, zorder=1)

                ax_i.set_xlim([0,1])
                ax_i.set_ylim([0,1])
                ax_i.set_xlabel('Allele frequency', fontsize = 14)
                ax_i.set_ylabel(r'$\mathrm{Pr}\left [ \mathrm{survival} \right ]$', fontsize = 14)

                ax_i.spines['top'].set_visible(False)

                line, = ax_i.plot([0,1],[1,1],'k:',linewidth=5)
                line.set_dashes((0.5, 0.5))

                ax_i.set_title( latex_dict[taxon], fontweight='bold', fontsize=17)

                #fig.suptitle(latex_dict[strain], fontsize=28, fontweight='bold')

            if (taxon_idx == 0):
                legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= 'Quasi-neutral'),
                                   Line2D([0], [0], ls=':', color='k', lw=1.5, label= 'Hitchhiking' )]
                ax_i.legend(handles=legend_elements, loc='upper left', fontsize=12)

    fig_name = pt.get_path() + '/figs/freq_vs_prob_fixation.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_allele_corr_delta(min_trajectory_length=3):

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    #transit_times = {}
    taxa = ['B', 'S', 'D']

    r2s_obs_dict = {}

    #r2s_null_dict = {}

    for treatment in ['0', '1']:
        for taxon in taxa:
            r2s = []
            for replicate in replicates:

                population = treatment + taxon + replicate
                if population in pt.populations_to_ignore:
                    continue
                sys.stderr.write("Processing %s...\n" % population)

                mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
                population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
                state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

                times = mutations[0][9]
                Ms = np.zeros_like(times)*1.0
                fixed_Ms = np.zeros_like(times)*1.0

                for mutation_idx_i in range(0,len(mutations)):

                    location_i, gene_name_i, allele_i, var_type_i, test_statistic_i, pvalue_i, cutoff_idx_i, depth_fold_change_i, depth_change_pvalue_i, times_i, alts_i, depths_i, clone_times_i, clone_alts_i, clone_depths_i = mutations[mutation_idx_i]

                    state_Ls_i = state_trajectories[mutation_idx_i]
                    good_idx_i, filtered_alts_i, filtered_depths_i = timecourse_utils.mask_timepoints(times_i, alts_i, depths_i, var_type_i, cutoff_idx_i, depth_fold_change_i, depth_change_pvalue_i)
                    freqs_i = timecourse_utils.estimate_frequencies(filtered_alts_i, filtered_depths_i)

                    masked_times_i = times[good_idx_i]
                    masked_freqs_i = freqs_i[good_idx_i]
                    masked_state_Ls_i = state_Ls_i[good_idx_i]

                    P_idx_i = np.where(masked_state_Ls_i == 3)[0]
                    if len(P_idx_i) < min_trajectory_length:
                        continue
                    first_P_i = P_idx_i[0]
                    last_P_i = P_idx_i[-1]

                    masked_freqs_P_i = masked_freqs_i[first_P_i:last_P_i+1]
                    masked_times_P_i = masked_times_i[first_P_i:last_P_i+1]

                    delta_masked_freqs_P_i = masked_freqs_P_i[1:] - masked_freqs_P_i[:-1]
                    delta_masked_times_P_i = masked_times_P_i[:-1]


                    for mutation_idx_j in range(mutation_idx_i+1,len(mutations)):

                        location_j, gene_name_j, allele_j, var_type_j, test_statistic_j, pvalue_j, cutoff_jdx_j, depth_fold_change_j, depth_change_pvalue_j, times_j, alts_j, depths_j, clone_times_j, clone_alts_j, clone_depths_j = mutations[mutation_idx_j]

                        state_Ls_j = state_trajectories[mutation_idx_j]
                        good_idx_j, filtered_alts_j, filtered_depths_j = timecourse_utils.mask_timepoints(times_j, alts_j, depths_j, var_type_j, cutoff_jdx_j, depth_fold_change_j, depth_change_pvalue_j)
                        freqs_j = timecourse_utils.estimate_frequencies(filtered_alts_j, filtered_depths_j)

                        masked_times_j = times[good_idx_j]
                        masked_freqs_j = freqs_j[good_idx_j]
                        masked_state_Ls_j = state_Ls_j[good_idx_j]

                        P_jdx_j = np.where(masked_state_Ls_j == 3)[0]
                        if len(P_jdx_j) < min_trajectory_length:
                          continue
                        first_P_j = P_jdx_j[0]
                        last_P_j = P_jdx_j[-1]

                        masked_freqs_P_j = masked_freqs_j[first_P_j:last_P_j+1]
                        masked_times_P_j = masked_times_j[first_P_j:last_P_j+1]

                        delta_masked_freqs_P_j = masked_freqs_P_j[1:] - masked_freqs_P_j[:-1]
                        # delta_f = f_t_plus_1 - f_t
                        delta_masked_times_P_j = masked_times_P_j[:-1]

                        intersect_times = np.intersect1d(delta_masked_times_P_i, delta_masked_times_P_j)

                        if len(intersect_times)>=3:

                            intersect_idx_i = [np.where(delta_masked_times_P_i == intersect_time)[0][0] for intersect_time in intersect_times ]
                            intersect_delta_i = delta_masked_freqs_P_i[intersect_idx_i]

                            intersect_idx_j = [np.where(delta_masked_times_P_j == intersect_time)[0][0] for intersect_time in intersect_times ]
                            intersect_delta_j = delta_masked_freqs_P_j[intersect_idx_j]



                            if len(intersect_delta_i) != len(intersect_delta_j):
                                print(len(intersect_delta_j), len(intersect_delta_j))

                            r2 = stats.pearsonr(intersect_delta_i, intersect_delta_j)[0] ** 2
                            r2s.append(r2)

            r2s_obs_dict[treatment + taxon] = r2s

    fig = plt.figure(figsize = (12, 12))

    tuples = [ (0,0), (0,1), (1,0)]
    for i, taxon in enumerate(taxa):
            ax_i = plt.subplot2grid((2, 2), tuples[i], colspan=1)

            for treatment in ['0', '1']:
                r2_treatment_taxon = r2s_obs_dict[treatment+taxon]

                ax_i.hist(r2_treatment_taxon, label=pt.get_treatment_name(treatment), color= pt.get_colors(treatment), lw=5, histtype='step', bins = 30, alpha =1, weights=np.zeros_like(r2_treatment_taxon) + 1. / len(r2_treatment_taxon))
                ax_i.set_xlim([0,1] )
                ax_i.set_yscale('log', basey=10)
                ax_i.set_title(latex_dict[taxon], fontsize=18, fontweight='bold')
                ax_i.set_xlabel("Squared correlation between\nallele frequency trajectories, " + r'$\rho_{M_{\mathrm{P}}^{ (i) }, M_{\mathrm{P}}^{ (j) }  }^{2} $' , fontsize = 14)
                ax_i.set_ylabel('Frequency', fontsize = 20 )

                if i == 0:
                    ax_i.legend(loc='upper right', fontsize=12)


    fig_name = pt.get_path() + '/figs/r2_B_S.png'
    fig.subplots_adjust(hspace=0.3)
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_0B3_big(population='0B3'):

    fig = plt.figure(figsize = (12, 6))

    annotated_timecourse_path = pt.get_path() + "/data/timecourse_final/%s_annotated_timecourse.txt" % population
    annotated_timecourse_file = open(annotated_timecourse_path ,"r")

    first_line = annotated_timecourse_file.readline()
    first_line = first_line.strip()
    first_line_items = first_line.split(",")
    times = np.asarray([float(x.strip().split(':')[1]) for x in first_line_items[13::2]])
    times = np.insert(times, 0, 0, axis=0)

    for i, line in enumerate(annotated_timecourse_file):
        line = line.strip()
        items = line.split(",")
        pass_or_fail = items[12].strip()
        if pass_or_fail == 'FAIL':
            continue

        alt_cov = np.asarray([ float(x) for x in items[13::2]])
        total_cov = np.asarray([ float(x) for x in items[14::2]])
        # pseudocount to avoid divide by zero error
        alt_cov = alt_cov
        total_cov = total_cov + 1
        freqs = alt_cov / total_cov
        freqs = np.insert(freqs, 0, 0, axis=0)

        rgb = pt.mut_freq_colormap()
        rgb = pt.lighten_color(rgb, amount=0.5)


        if len(times) == len(freqs) + 1:
            freqs = np.insert(freqs, 0, 0, axis=0)

        plt.plot(times, freqs, '.-', c=rgb, lw=3, alpha=0.4)


    plt.xlim([0,max(times)])
    plt.ylim([0, 1])

    plt.ylabel( str(10**int(population[0])) + '-day transfers', fontsize =12  )

    plt.tick_params(axis="x", labelsize=8)
    plt.tick_params(axis="y", labelsize=8)


    fig.text(0.5, 0.04, 'Days, ' + r'$t$', ha='center', va='center', fontsize=24)
    fig.text(0.05, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=24)


    fig.suptitle(latex_dict['B'], fontsize=28, fontweight='bold')
    fig_name = pt.get_path() + '/figs/mut_trajectories_0B3.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_B_S_multiplicity():

    taxa = ['B', 'S']

    parallelism_axes = {}

    fig = plt.figure(figsize = (12, 12))

    gene_data = parse_file.parse_gene_list('B')

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    # to get the common gene names for each ID

    for treatment_idx, treatment in enumerate(treatments):

        ax_multiplicity = plt.subplot2grid((3, 3), (0, treatment_idx), colspan=1)
        ax_regression = plt.subplot2grid((3, 3), (1, treatment_idx), colspan=1)
        ax_venn = plt.subplot2grid((3, 3), (2, treatment_idx), colspan=1)

        ax_multiplicity.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
        ax_multiplicity.set_ylabel('Fraction mutations ' + r'$\geq m$', fontsize=14)
        ax_multiplicity.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

        ax_multiplicity.set_xscale('log', basex=10)
        ax_multiplicity.set_yscale('log', basey=10)

        ax_multiplicity.set_ylim([0.001, 1.1])
        ax_multiplicity.set_xlim([0.07, 130])

        ax_regression.set_xlabel('Gene multiplicity, ' + r'$m$' + '\n' + latex_dict['B']   , fontsize=14)
        ax_regression.set_ylabel('Gene multiplicity, ' + r'$m$' + '\n' + latex_dict['S'] , fontsize=14)

        ax_regression.set_xscale('log', basex=10)
        ax_regression.set_yscale('log', basey=10)

        ax_venn.axis('off')

        mult_taxa_dict = {}

        for taxon in taxa:

            if taxon == 'B':
                multiplicity_label = r'$\mathrm{wt}$'
            else:
                multiplicity_label = r'$\Delta \mathrm{spo0A}$'

            populations = [treatment+taxon + replicate for replicate in replicates ]

            # Load convergence matrix
            convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

            gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

            G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

            sys.stdout.write("Total parallelism for %s = %g (p=%g)\n" % (treatment+taxon, G,pvalue))

            predictors = []
            responses = []

            gene_hits = []
            gene_predictors = []

            Ls = []

            for gene_name in convergence_matrix.keys():

                # get multiplicities for regression
                if gene_name not in mult_taxa_dict:
                    mult_taxa_dict[gene_name] = {}
                    mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']
                else:
                    mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']

                convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

                Ls.append(convergence_matrix[gene_name]['length'])
                m = gene_parallelism_statistics[gene_name]['multiplicity']

                n = 0
                nfixed = 0

                for population in populations:
                    for t,L,f in convergence_matrix[gene_name]['mutations'][population]:
                        fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                        predictors.append(m)
                        responses.append(fixed_weight)

                        n+=1
                        nfixed+=fixed_weight

                if n > 0.5:
                    gene_hits.append(n)
                    gene_predictors.append(m)

            Ls = np.asarray(Ls)
            ntot = len(predictors)
            mavg = ntot*1.0/len(Ls)

            predictors, responses = (np.array(x) for x in zip(*sorted(zip(predictors, responses), key=lambda pair: (pair[0]))))

            gene_hits, gene_predictors = (np.array(x) for x in zip(*sorted(zip(gene_hits, gene_predictors), key=lambda pair: (pair[0]))))

            rescaled_predictors = np.exp(np.fabs(np.log(predictors/mavg)))

            #logit_mod = sm.Logit(responses, sm.add_constant(rescaled_predictors))
            #logit_res = logit_mod.fit()

            #sys.stdout.write("Logistic regression for %ss:\n" % (population))
            #sys.stdout.write("%s\n" % str(logit_res.summary()))
            #sys.stdout.write("Params:\n")
            #sys.stdout.write("%s\n" % str(logit_res.params))

            sys.stdout.write("Avg mutation multiplicity=%g, Avg fixed mutation multiplicity=%g\n" % (predictors.sum()/len(responses), (predictors*responses).sum()/responses.sum()))
            sys.stderr.write("Calculating null distribution...\n")
            null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)

            # default base is 10
            theory_ms = np.logspace(-2,2,100)
            theory_survivals = null_survival_function(theory_ms)
            theory_survivals /= theory_survivals[0]

            sys.stderr.write("Done!\n")

            # step function
            ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3, color=pt.get_colors(treatment),alpha=0.8, ls=pt.get_taxon_ls(taxon), label='Observed ' + multiplicity_label, drawstyle='steps', zorder=2)

            ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color='grey',alpha=0.8, ls=pt.get_taxon_ls(taxon), label= 'Null ' +  multiplicity_label, zorder=1)

        if treatment_idx == 0:
            ax_multiplicity.legend( loc='lower left', fontsize=8)

        for gene_name, gene_dict in mult_taxa_dict.items():
            if 'B' not in gene_dict:
                mult_taxa_dict[gene_name]['B'] = 0
            if 'S' not in gene_dict:
                mult_taxa_dict[gene_name]['S'] = 0

        #mult_B_S = [(mult_taxa_dict[gene_name]['B'], mult_taxa_dict[gene_name]['S']) for gene_name in sorted(mult_taxa_dict) if (mult_taxa_dict[gene_name]['B'] > 0) and (mult_taxa_dict[gene_name]['S'] > 0) ]

        # then get venn diagram
        # import significant genes
        parallel_genes_B = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+'B'), "r")
        parallel_genes_S = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+'S'), "r")

        gene_significant_multiplicity_dict = {}


        parallel_genes_B_list = []
        for i, line in enumerate(parallel_genes_B):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            parallel_genes_B_list.append(items[0])

            if items[0] not in gene_significant_multiplicity_dict:
                gene_significant_multiplicity_dict[items[0]] = {}
            gene_significant_multiplicity_dict[items[0]]['B'] = float(items[6].strip())

        parallel_genes_S_list = []
        for i, line in enumerate(parallel_genes_S):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            parallel_genes_S_list.append(items[0])

            if items[0] not in gene_significant_multiplicity_dict:
                gene_significant_multiplicity_dict[items[0]] = {}
            gene_significant_multiplicity_dict[items[0]]['S'] = float(items[6].strip())


        #mult_B_S = [(mult_taxa_dict[gene_name]['B'], mult_taxa_dict[gene_name]['S']) for gene_name in sorted(mult_taxa_dict) if (mult_taxa_dict[gene_name]['B'] > 0) and (mult_taxa_dict[gene_name]['S'] > 0) ]
        mult_B_S = [(gene_significant_multiplicity_dict[gene_name]['B'], gene_significant_multiplicity_dict[gene_name]['S']) for gene_name in sorted(gene_significant_multiplicity_dict) if ('B' in gene_significant_multiplicity_dict[gene_name]) and ('S' in gene_significant_multiplicity_dict[gene_name]) ]

        mult_B = [x[0] for x in mult_B_S]
        mult_S = [x[1] for x in mult_B_S]


        ax_regression.scatter(mult_B, mult_S, color=pt.get_colors(treatment),alpha=1,s=90, zorder=2)

        if treatment_idx == 2:
            ax_regression.set_xlim([  0.9, max( mult_B + mult_S )*1.5  ])
            ax_regression.set_ylim([  0.9, max( mult_B + mult_S )*1.5  ])
            ax_regression.plot([0.9, max( mult_B + mult_S )*1.5 ], [ 0.9, max( mult_B + mult_S )*1.5  ], lw = 3, c='k', ls = '--', zorder=1 )


        else:
            ax_regression.set_xlim([  min( mult_B + mult_S )/1.5, max( mult_B + mult_S )*1.5  ])
            ax_regression.set_ylim([  min( mult_B + mult_S )/1.5, max( mult_B + mult_S )*1.5  ])
            ax_regression.plot([  min( mult_B + mult_S )/1.5, max( mult_B + mult_S )*1.5 ], [  min( mult_B + mult_S )/1.5, max( mult_B + mult_S )*1.5  ], lw = 3, c='k', ls = '--', zorder=1 )

        venn = venn2(subsets = (len(parallel_genes_B_list), len(parallel_genes_S_list), len(set(parallel_genes_B_list) & set(parallel_genes_S_list))), ax=ax_venn, set_labels=('', ''), set_colors=(pt.get_colors(treatment), pt.get_colors(treatment)))
        c = venn2_circles(subsets=(len(parallel_genes_B_list), len(parallel_genes_S_list), len(set(parallel_genes_B_list) & set(parallel_genes_S_list))), ax=ax_venn, linestyle='dashed')
        #set_colors=(pt.get_colors(treatment), pt.get_colors(treatment)),

        c[0].set_ls('--')
        c[1].set_ls(':')
        c[0].set_lw(5)
        c[1].set_lw(5)
        c[0].set_edgecolor(pt.get_colors(treatment))
        c[1].set_edgecolor(pt.get_colors(treatment))
        #c[0].set_radius(len(parallel_genes_B_list) / 80 )
        #c[1].set_radius(len(parallel_genes_S_list) / 80)

    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    fig_name = pt.get_path() + '/figs/B_S_multiplicity.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_within_taxon_paralleliism(taxon):

    fig = plt.figure(figsize = (18, 12))

    gene_data = parse_file.parse_gene_list('B')

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    # to get the common gene names for each ID

    ax_multiplicity = plt.subplot2grid((2, 3), (0, 0), colspan=1)
    ax_venn = plt.subplot2grid((2, 3), (0, 1), colspan=1)
    #ax_mult_pfix = plt.subplot2grid((3, 2), (1, 0), colspan=1)
    ax_mult_freq = plt.subplot2grid((2, 3), (0, 2), colspan=1)
    ax_mult_1_10 = plt.subplot2grid((2, 3), (1, 0), colspan=1)

    ax

    ax_multiplicity.set_xscale('log', basex=10)
    ax_multiplicity.set_yscale('log', basey=10)
    ax_multiplicity.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_multiplicity.set_ylabel('Fraction mutations ' + r'$\geq m$', fontsize=14)

    ax_multiplicity.set_ylim([0.001, 1.1])
    ax_multiplicity.set_xlim([0.07, 130])

    ax_venn.axis('off')

    ax_mult_pfix.set_xscale('log', basex=10)
    ax_mult_pfix.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)

    ax_mult_freq.set_xscale('log', basex=10)
    ax_mult_freq.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)

    significant_multiplicity_dict = {}

    for treatment_idx, treatment in enumerate(treatments):

        significan_multiplicity_taxon = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon), "r")
        significan_multiplicity_list = []
        for i, line in enumerate(significan_multiplicity_taxon):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            significan_multiplicity_list.append(items[0])

        significant_multiplicity_dict[treatment] = significan_multiplicity_list

        populations = [treatment+taxon + replicate for replicate in replicates ]

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        sys.stdout.write("Total parallelism for %s = %g (p=%g)\n" % (treatment+taxon, G,pvalue))

        predictors = []
        responses = []

        gene_hits = []
        gene_predictors = []

        Ls = []

        for gene_name in convergence_matrix.keys():

            # try to get probability of fixation for all treatments with at least 30 total mutations

            # get multiplicities for regression
            #if gene_name not in mult_taxa_dict:
            #    mult_taxa_dict[gene_name] = {}
            #    mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']
            #else:
            #    mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']

            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

            Ls.append(convergence_matrix[gene_name]['length'])
            m = gene_parallelism_statistics[gene_name]['multiplicity']

            n = 0
            nfixed = 0

            for population in populations:
                for t,L,f in convergence_matrix[gene_name]['mutations'][population]:
                    fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                    predictors.append(m)
                    responses.append(fixed_weight)

                    n+=1
                    nfixed+=fixed_weight

            if n > 0.5:
                gene_hits.append(n)
                gene_predictors.append(m)

        Ls = np.asarray(Ls)
        ntot = len(predictors)
        mavg = ntot*1.0/len(Ls)

        predictors, responses = (np.array(x) for x in zip(*sorted(zip(predictors, responses), key=lambda pair: (pair[0]))))

        gene_hits, gene_predictors = (np.array(x) for x in zip(*sorted(zip(gene_hits, gene_predictors), key=lambda pair: (pair[0]))))

        rescaled_predictors = np.exp(np.fabs(np.log(predictors/mavg)))

        null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)

        # default base is 10
        theory_ms = np.logspace(-2,2,100)
        theory_survivals = null_survival_function(theory_ms)
        theory_survivals /= theory_survivals[0]

        sys.stderr.write("Done!\n")

        # step function
        ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3.5, color=pt.get_colors(treatment),alpha=0.8, ls=pt.get_taxon_ls(taxon), label= str(int(10**int(treatment))) + '-day', drawstyle='steps', zorder=2)

        ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color=pt.get_colors(treatment), alpha=0.8, ls=':',  zorder=1)

    print(len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])))

    subset_tuple = (len( significant_multiplicity_dict['0']), \
                    len( significant_multiplicity_dict['1']), \
                    len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['1'])), \
                    len(significant_multiplicity_dict['2']), \
                    len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['2'])), \
                    len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])),  \
                    len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])))


    venn = venn3(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1'), pt.get_colors('2')))
    c = venn3_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')

    fig.suptitle(latex_dict[taxon], fontsize=30)

    fig.subplots_adjust() #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/multiplicity_%s.png" % taxon
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






#plot_spo0a_fitness()



#plot_spores()

#plot_within_taxon_paralleliism('D')

#plot_B_S_multiplicity()

#plot_0B3_big()
#plot_allele_corr_delta()

#allele_survival()
#plot_B_S_mutation_trajectory()

#get_mult_similarity()

#plot_allele_freqs_all_treats('B')
#plot_allele_freqs_all_treats('S')
#plot_allele_freqs_all_treats('C')

#plot_taxon_mutation_trajectory('C')


#plot_bPTR_all()
temporal_coverage_B_S()
