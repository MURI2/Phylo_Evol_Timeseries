from __future__ import division
import os, sys, json
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.colors import ColorConverter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


import phylo_tools as pt
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.stats.api as sms
import statsmodels.formula.api as smf
from statsmodels.compat import lzip

#from sklearn.metrics import pairwise_distances
#from skbio.stats.ordination import pcoa

import parse_file
import timecourse_utils
import mutation_spectrum_utils

from scipy.special import gammaln

from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

from sklearn.decomposition import PCA


np.random.seed(123456789)

treatments=pt.treatments
replicates = pt.replicates


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']


#legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}_{WT} (t)$'),
#                   Line2D([0], [0], ls=':', color='k', lw=1.5, label= r'$\overline{M}_{\Delta \mathrm{spo0A}} (t)$')]
#ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=12)


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

    fig_name = pt.get_path() + '/figs/spore_assay.pdf'
    fig.savefig(fig_name, format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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

    # mean initia CFU counts from notebook

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

    fig_name = pt.get_path() + '/figs/fitness_spo0a.pdf'
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def temporal_plasmid_coverage_B_S():

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

                times_to_ignore = pt.samples_to_remove(population)

                for time in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:

                    sample = '%s_%s' % (population, time)
                    sample_path = pt.get_path() + '/data/rebreseq_json/%s.json' % sample

                    if os.path.isfile(sample_path) == False:
                        continue

                    if times_to_ignore != None:
                        if time in times_to_ignore:
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
    fig.text(0.04, 0.5, 'Plasmid-chromosome coverage ratio', va='center', rotation='vertical', fontsize=28)

    fig_name = pt.get_path() + '/figs/plasmid_coverage_B_S.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






def temporal_plasmid_coverage_D():

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

                    #if population in pt.samples_to_remove(population)

                    times_to_ignore = pt.samples_to_remove(population)
                    if times_to_ignore != None:
                        if time in times_to_ignore:
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
    fig.text(0.05, 0.5, 'Plasmid-chromosome 1 coverage ratio', va='center', rotation='vertical', fontsize=28)

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

            population = treatment + strain + replicate
            annotated_timecourse_path = pt.get_path() + "/data/timecourse_final/%s_annotated_timecourse.txt" % population
            if os.path.exists(annotated_timecourse_path) == False:
                continue
            annotated_timecourse_file = open(annotated_timecourse_path ,"r")
            print(population)


            ax_i = plt.subplot2grid((3, 5), (row_count, column_count), colspan=1)

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
    #fig_name = pt.get_path() + '/figs/mut_trajectories_%s.png'
    fig_name = pt.get_path() + '/figs/mut_trajectories_%s.pdf'
    #fig.savefig(fig_name % strain, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    fig.savefig(fig_name % strain, bbox_inches = "tight",format='pdf', pad_inches = 0.4, dpi = 600)

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





def plot_mutation_trajectory_B_S():
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

    # fixed mutations t-test for B vs S at day 1 and 10 at day 600

    for treatment in ['0', '1']:
        fixation_count_dict = {}
        for taxon in taxa:
            fixation_count_dict[taxon] = []
            for replicate in replicates:
                population = treatment + taxon + replicate
                fixation_idx = np.where(fixed_mutation_trajectories[population][0] == 600)
                #print(fixed_mutation_trajectories[population][1])
                fixation_count_dict[taxon].append(fixed_mutation_trajectories[population][1][fixation_idx][0])
        #print(fixation_count_dict)

        #print(np.mean(fixation_count_dict['B']))

        #print(stats.ttest_ind(fixation_count_dict['B'], fixation_count_dict['S'], equal_var=False))

        #ss.t.sf(t_stat, len(x) + len(y) - 2)

    fig = plt.figure(figsize = (10, 9))

    column_count = 0

    for treatment in treatments:

        ax_t_vs_M = plt.subplot2grid((3, 3), (0, column_count), colspan=1)

        ax_t_vs_delta_M = plt.subplot2grid((3, 3), (1, column_count), colspan=1)

        ax_M_vs_F = plt.subplot2grid((3, 3), (2, column_count), colspan=1)

        ax_t_vs_M.text(-0.1, 1.07, sub_plot_labels[column_count], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_t_vs_M.transAxes)
        ax_t_vs_delta_M.text(-0.1, 1.07, sub_plot_labels[column_count+3], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_t_vs_delta_M.transAxes)
        ax_M_vs_F.text(-0.1, 1.07, sub_plot_labels[column_count+6], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_M_vs_F.transAxes)


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
                ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=12)

        ax_t_vs_M.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

        if treatment == '2':
            ax_M_vs_F.yaxis.set_major_locator(MaxNLocator(integer=True))

        if column_count == 0:

            ax_t_vs_M.set_ylabel('Mutations, ' + r'$M(t)$', fontsize = 15)
            ax_M_vs_F.set_ylabel('Fixed mutations', fontsize = 15)
            ax_t_vs_delta_M.set_ylabel('Change in mutations,\n' + r'$M(t)/M(t-1)$', fontsize = 15)

        ax_t_vs_M.set_ylim([0.008 , 300])
        ax_t_vs_delta_M.set_ylim([0.1 , 70])
        #ax_M_vs_F.set_ylim([-0.1 , 44])

        column_count += 1

    fig.text(0.53, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=28)


    fig_name = pt.get_path() + '/figs/rate_B_S.pdf'
    fig.savefig(fig_name, format='pdf',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_mutation_trajectory_taxon(taxon):

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

        ax_t_vs_M.text(-0.1, 1.07, sub_plot_labels[column_count], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_t_vs_M.transAxes)
        ax_t_vs_delta_M.text(-0.1, 1.07, sub_plot_labels[column_count+3], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_t_vs_delta_M.transAxes)
        ax_M_vs_F.text(-0.1, 1.07, sub_plot_labels[column_count+6], fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_M_vs_F.transAxes)


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


        # keep them on the same y axes
        if taxon == 'C':
            ax_t_vs_delta_M.set_ylim([0.2,42])
        elif taxon == 'D':
            ax_t_vs_delta_M.set_ylim([0.2,20])



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


    fig_name = pt.get_path() + '/figs/rate_%s.pdf' % taxon
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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



def plot_allele_corr_delta_B_S(min_trajectory_length=3):

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    #transit_times = {}
    taxa = ['B', 'S']

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

    fig = plt.figure(figsize = (12, 6))

    #tuples = [ (0,0), (0,1), (1,0), (1,1)]
    #for i, taxon in enumerate(taxa):
    for treatment_idx, treatment in enumerate(['0', '1']):

        ax_i = plt.subplot2grid((1, 2), (0,treatment_idx), colspan=1)
        ax_i.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)
        ax_i.text(-0.1, 1.07, sub_plot_labels[treatment_idx], fontsize=22, fontweight='bold', ha='center', va='center', transform=ax_i.transAxes)

        print(treatment)
        print(stats.ks_2samp(r2s_obs_dict[treatment+'B'], r2s_obs_dict[treatment+'S']))

        for taxon in taxa:

            r2_treatment_taxon = r2s_obs_dict[treatment+taxon]
            ax_i.hist(r2_treatment_taxon, label=latex_dict[taxon], linestyle=pt.get_taxon_ls(taxon), color= pt.get_colors(treatment), lw=3, histtype='step', bins = 40, alpha =1, weights=np.zeros_like(r2_treatment_taxon) + 1. / len(r2_treatment_taxon))
            ax_i.set_xlim([0,1] )
            ax_i.set_yscale('log', basey=10)
            #ax_i.set_title(latex_dict[taxon], fontsize=18, fontweight='bold')
            #ax_i.set_xlabel("Squared correlation between\nallele frequency trajectories, " + r'$\rho_{M_{\mathrm{P}}^{ (i) }, M_{\mathrm{P}}^{ (j) }  }^{2} $' , fontsize = 14)

        if treatment_idx == 0:
            ax_i.legend(loc='upper right', fontsize=12)


    fig.text(0.5, -0.03, "Squared correlation between\nallele frequency trajectories, " + r'$\rho_{M_{\mathrm{P}}^{ (i) }, M_{\mathrm{P}}^{ (j) }  }^{2} $', ha='center', va='center', fontsize=20)
    fig.text(0.05, 0.5, 'Frequency', ha='center', va='center', rotation='vertical',  fontsize=24)

    fig_name = pt.get_path() + '/figs/r2_B_S.pdf'
    fig.subplots_adjust(hspace=0.3)
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_allele_corr_delta(min_trajectory_length=3):

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    #transit_times = {}
    taxa = ['C', 'D', 'J', 'P', 'F']

    r2s_obs_dict = {}

    for taxon in taxa:
        for treatment in ['0', '1']:
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

    fig = plt.figure(figsize = (12, 6))

    tuples = [ (0,0), (0,1), (0,2), (1,0), (1,1)]
    #for treatment_idx, treatment in enumerate(['0', '1']):

    sub_plot_count = 0

    for taxon_idx, taxon in enumerate(taxa):
        ax_i = plt.subplot2grid((2, 3), tuples[taxon_idx], colspan=1)
        ax_i.set_title(latex_dict[taxon], fontsize=13)
        ax_i.text(-0.1, 1.07, sub_plot_labels[sub_plot_count], fontsize=20, fontweight='bold', ha='center', va='center', transform=ax_i.transAxes)

        sub_plot_count+=1

        for treatment in ['0', '1']:
            r2_treatment_taxon = r2s_obs_dict[treatment+taxon]
            if len(r2_treatment_taxon) == 0:
                print("No correlations!")
                continue
            ax_i.hist(r2_treatment_taxon, label=str(10**int(population[0])) + '-day transfers', linestyle=pt.get_taxon_ls(taxon), color= pt.get_colors(treatment), lw=3, histtype='step', bins = 40, alpha=1, weights=np.zeros_like(r2_treatment_taxon) + 1. / len(r2_treatment_taxon))
            ax_i.set_xlim([0,1] )
            ax_i.set_ylim([0.01,0.2] )
            ax_i.set_yscale('log', basey=10)


        if taxon_idx == 0:
            ax_i.legend(loc='upper right', fontsize=12)

    fig.text(0.5, 0.03, "Squared correlation between allele frequency trajectories, " + r'$\rho_{M_{\mathrm{P}}^{ (i) }, M_{\mathrm{P}}^{ (j) }  }^{2} $', ha='center', va='center', fontsize=20)
    fig.text(0.05, 0.5, 'Frequency', ha='center', va='center', rotation='vertical',  fontsize=24)

    fig_name = pt.get_path() + '/figs/r2_all.pdf'
    fig.subplots_adjust(hspace=0.3)
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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



def plot_B_S_paralleliism():

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
        ax_multiplicity.text(-0.1, 1.07, sub_plot_labels[treatment_idx], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)


        ax_regression.set_xlabel('Gene multiplicity, ' + r'$m$' + '\n' + latex_dict['B']   , fontsize=14)
        ax_regression.set_ylabel('Gene multiplicity, ' + r'$m$' + '\n' + latex_dict['S'] , fontsize=14)

        ax_regression.set_xscale('log', basex=10)
        ax_regression.set_yscale('log', basey=10)
        ax_regression.text(-0.1, 1.07, sub_plot_labels[3+treatment_idx], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_regression.transAxes)


        ax_venn.axis('off')
        ax_venn.text(-0.1, 1.07, sub_plot_labels[6+treatment_idx], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_venn.transAxes)


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
                    for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:
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

            #rescaled_predictors = np.exp(np.fabs(np.log(predictors/mavg)))

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


        mult_B_S = [(mult_taxa_dict[gene_name]['B'], mult_taxa_dict[gene_name]['S']) for gene_name in sorted(mult_taxa_dict) if (mult_taxa_dict[gene_name]['B'] > 0) and (mult_taxa_dict[gene_name]['S'] > 0) ]
        mult_significant_B_S = [(gene_significant_multiplicity_dict[gene_name]['B'], gene_significant_multiplicity_dict[gene_name]['S']) for gene_name in sorted(gene_significant_multiplicity_dict) if ('B' in gene_significant_multiplicity_dict[gene_name]) and ('S' in gene_significant_multiplicity_dict[gene_name]) ]

        mult_B = [x[0] for x in mult_B_S]
        mult_S = [x[1] for x in mult_B_S]

        mult_significant_B = [x[0] for x in mult_significant_B_S]
        mult_significant_S = [x[1] for x in mult_significant_B_S]


        ax_regression.scatter(mult_B, mult_S, color=pt.get_colors(treatment),alpha=1,s=90, zorder=2)

        ax_regression.scatter(mult_significant_B, mult_significant_S, linewidth=3, facecolors=pt.get_colors(treatment), edgecolors='k', alpha=1, s=90, zorder=3)

        ax_regression.set_xlim([  0.05, 200 ])
        ax_regression.set_ylim([  0.05, 200 ])
        ax_regression.plot([   0.02, 200  ], [   0.02, 200   ], lw = 3, c='k', ls = '--', zorder=1 )


        if len(mult_significant_B) >= 3:
            slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mult_significant_B), np.log10(mult_significant_S))

            # hypothetical slope of 1
            ratio = (slope - 1) / std_err
            pval = stats.t.sf(np.abs(ratio), len(mult_significant_B)-2)*2
            # two sided or one sided?
            sys.stderr.write("Treatment %s-day slope t-test = %f, p = %f, df=%d\n" %  (str(10**int(treatment)), round(ratio, 3),  round(pval, 3),  len(mult_significant_B)-2  ) )


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
    fig_name = pt.get_path() + '/figs/multiplicity_B_vs_S.pdf'
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_within_taxon_paralleliism(taxon):

    fig = plt.figure(figsize = (12, 8))

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    # to get the common gene names for each ID

    ax_multiplicity = plt.subplot2grid((2, 3), (0, 0), colspan=1)
    ax_mult_freq = plt.subplot2grid((2, 3), (0, 1), colspan=1)
    ax_venn = plt.subplot2grid((2, 3), (0, 2), colspan=1)

    ax_multiplicity.set_xscale('log', basex=10)
    ax_multiplicity.set_yscale('log', basey=10)
    ax_multiplicity.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_multiplicity.set_ylabel('Fraction mutations ' + r'$\geq m$', fontsize=14)
    ax_multiplicity.text(-0.1, 1.07, sub_plot_labels[0], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)


    ax_multiplicity.set_ylim([0.001, 1.1])
    ax_multiplicity.set_xlim([0.07, 130])


    ax_mult_freq.set_xscale('log', basex=10)
    ax_mult_freq.set_yscale('log', basey=10)
    ax_mult_freq.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_mult_freq.set_ylabel('Mean maximum allele frequency, ' + r'$f_{max}$', fontsize=11)
    ax_mult_freq.text(-0.1, 1.07, sub_plot_labels[1], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_mult_freq.transAxes)


    ax_venn.axis('off')
    ax_venn.text(-0.1, 1.07, sub_plot_labels[2], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_venn.transAxes)


    alpha_treatment_dict = {'0':0.5, '1':0.5, '2':0.8}


    significant_multiplicity_dict = {}

    significant_multiplicity_values_dict = {}

    multiplicity_dict = {}

    g_score_p_label_dict = {}

    all_mults = []
    all_freqs = []

    treatments_in_taxon = []

    label_y_axes = [0.3,0.2,0.1]

    for treatment_idx, treatment in enumerate(treatments):


        significan_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significan_multiplicity_taxon_path) == False:
            continue
        treatments_in_taxon.append(treatment)
        significan_multiplicity_taxon = open(significan_multiplicity_taxon_path, "r")

        significan_multiplicity_list = []
        for i, line in enumerate( significan_multiplicity_taxon ):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            significan_multiplicity_list.append(items[0])

            if items[0] not in significant_multiplicity_values_dict:
                significant_multiplicity_values_dict[items[0]] = {}
                significant_multiplicity_values_dict[items[0]][treatment] = float(items[-2])
            else:
                significant_multiplicity_values_dict[items[0]][treatment] = float(items[-2])


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
        mean_gene_freqs = []

        Ls = []

        ax_mult_freqs_x = []
        ax_mult_freqs_y = []

        for gene_name in convergence_matrix.keys():

            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

            Ls.append(convergence_matrix[gene_name]['length'])
            m = gene_parallelism_statistics[gene_name]['multiplicity']


            if gene_name not in multiplicity_dict:
                multiplicity_dict[gene_name] = {}
                multiplicity_dict[gene_name][treatment] = m
            else:
                multiplicity_dict[gene_name][treatment] = m

            n = 0
            nfixed = 0
            freqs = []
            nf_max = 0

            for population in populations:
                for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:
                    fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                    predictors.append(m)
                    responses.append(fixed_weight)

                    n+=1
                    nfixed+=fixed_weight

                    # get freqs for regression
                    #if L == parse_file.POLYMORPHIC:
                    #freqs.append(f_max)
                    nf_max+=timecourse_utils.calculate_fixed_weight(L,f_max)

            if n > 0.5:
                gene_hits.append(n)
                gene_predictors.append(m)
                #mean_gene_freqs.append(np.mean(freqs))

                if nf_max > 0:
                    ax_mult_freqs_x.append(m)
                    ax_mult_freqs_y.append( nf_max / n )

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

        ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color=pt.get_colors(treatment), alpha=0.8, ls=':',  zorder=1)

        ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3, color=pt.get_colors(treatment),alpha=0.8, ls='--', label= str(int(10**int(treatment))) + '-day', drawstyle='steps', zorder=2)

        #ax_multiplicity.text(0.2, 0.3, g_score_p_label_dict['0'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)
        #ax_multiplicity.text(0.2, 0.2, g_score_p_label_dict['1'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)
        #ax_multiplicity.text(0.2, 0.1, g_score_p_label_dict['2'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)

        if pvalue < 0.001:
            pretty_pvalue = r'$\ll 0.001$'
        else:
            pretty_pvalue = '=' + str(round( pvalue, 4 ))

        g_score_p_label = r'$\Delta \ell_{{{}}}=$'.format(str(10**int(treatment))) + str(round(G, 3)) + ', ' + r'$P$' + pretty_pvalue

        text_color = pt.lighten_color(pt.get_colors(treatment), amount=1.3)

        ax_multiplicity.text(0.26, label_y_axes[treatment_idx], g_score_p_label, fontsize=7, ha='center', va='center',  color='k', transform=ax_multiplicity.transAxes)



        ax_mult_freq.scatter(ax_mult_freqs_x, ax_mult_freqs_y, color=pt.get_colors(treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=alpha_treatment_dict[treatment])

        all_mults.extend(ax_mult_freqs_x)
        all_freqs.extend(ax_mult_freqs_y)

        #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(ax_mult_freqs_x), np.log10(ax_mult_freqs_y))
        #print(slope, p_value)

    # make treatment pairs
    treatments_in_taxon.sort(key=float)

    for i in range(0, len(treatments_in_taxon)):

        for j in range(i+1, len(treatments_in_taxon)):

            ax_mult_i_j = plt.subplot2grid((2, 3), (1, i+j-1), colspan=1)
            ax_mult_i_j.set_xscale('log', basex=10)
            ax_mult_i_j.set_yscale('log', basey=10)
            ax_mult_i_j.set_xlabel( str( 10** int( treatments_in_taxon[i] ) ) + '-day gene multiplicity, ' + r'$m$', fontsize=14)
            ax_mult_i_j.set_ylabel(str( 10** int( treatments_in_taxon[j] ) ) + '-day gene multiplicity, ' + r'$m$', fontsize=14)
            ax_mult_i_j.plot([   0.05, 200  ], [   0.05, 200   ], lw = 3, c='k', ls = '--', zorder=1 )
            ax_mult_i_j.set_xlim([   0.05, 200  ])
            ax_mult_i_j.set_ylim([   0.05, 200  ])

            ax_mult_i_j.text(-0.1, 1.07, sub_plot_labels[2+i+j], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_mult_i_j.transAxes)


            multiplicity_pair = [(multiplicity_dict[gene_name][treatments_in_taxon[i]], multiplicity_dict[gene_name][treatments_in_taxon[j]]) for gene_name in sorted(multiplicity_dict) if (multiplicity_dict[gene_name][treatments_in_taxon[i]] > 0) and (multiplicity_dict[gene_name][treatments_in_taxon[j]] > 0) ]
            significant_multiplicity_pair = [(significant_multiplicity_values_dict[gene_name][treatments_in_taxon[i]], significant_multiplicity_values_dict[gene_name][treatments_in_taxon[j]]) for gene_name in sorted(significant_multiplicity_values_dict) if (treatments_in_taxon[i] in significant_multiplicity_values_dict[gene_name]) and (treatments_in_taxon[j] in significant_multiplicity_values_dict[gene_name]) ]

            # get mean colors
            ccv = ColorConverter()

            color_1 = np.array(ccv.to_rgb( pt.get_colors( treatments_in_taxon[i] ) ))
            color_2 = np.array(ccv.to_rgb( pt.get_colors( treatments_in_taxon[j] ) ))

            mix_color = 0.7 * (color_1 + color_2)
            mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)

            if (treatments_in_taxon[i] == '0') and (treatments_in_taxon[j] == '1'):
                #mix_color = pt.lighten_color(mix_color, amount=2.8)
                mix_color = 'gold'

            ax_mult_i_j.scatter([x[0] for x in multiplicity_pair], [x[1] for x in multiplicity_pair], marker=pt.plot_species_marker(taxon), facecolors=mix_color, edgecolors='none', alpha=0.8,s=90, zorder=2)
            ax_mult_i_j.scatter([x[0] for x in significant_multiplicity_pair], [x[1] for x in significant_multiplicity_pair], marker=pt.plot_species_marker(taxon), facecolors=mix_color, edgecolors='k', lw=1.5, alpha=1,s=90, zorder=3)

    if taxon == 'F':
        subset_tuple = (len( significant_multiplicity_dict['0']), \
                        len( significant_multiplicity_dict['1']), \
                        len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['1'])))

        venn = venn2(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1')))
        c = venn2_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')


    elif taxon == 'J':
        subset_tuple = (len( significant_multiplicity_dict['0']), \
                        len(significant_multiplicity_dict['2']), \
                        len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['2'])))

        venn = venn2(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('2')))
        c = venn2_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')


    else:
        subset_tuple = (len( significant_multiplicity_dict['0']), \
                        len( significant_multiplicity_dict['1']), \
                        len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['1'])), \
                        len(significant_multiplicity_dict['2']), \
                        len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['2'])), \
                        len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])),  \
                        len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])))

        venn = venn3(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1'), pt.get_colors('2')))
        c = venn3_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')


    ax_mult_freq.set_xlim([min(all_mults)*0.5, max(all_mults)*1.5])
    ax_mult_freq.set_ylim([min(all_freqs)*0.5, max(all_freqs)*1.5])




    fig.suptitle(latex_dict[taxon], fontsize=30)

    fig.subplots_adjust(wspace=0.3) #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/multiplicity_%s.pdf" % taxon
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()









def mutation_regression_B_S(time_measure = 'generations', set_time=500, slope_null=1):

    mutation_trajectories = {}

    taxa = ['B', 'S']
    # loop through taxa and get M(700) for all reps in each treatment
    for taxon in taxa:
        for treatment in treatments:
            for replicate in replicates:

                population = treatment + taxon + replicate
                if population in pt.populations_to_ignore:
                    continue

                sys.stderr.write("Processing %s...\t" % population)

                times, Ms, fixed_Ms = get_mutation_fixation_trajectories(population)

                if set_time not in times:
                    continue
                time_idx = np.where(times == set_time)

                #mutation_trajectories[population] = (times,np.log10(Ms))
                mutation_trajectories[population] = (times[time_idx][0],np.log10(Ms[time_idx][0] ))

                sys.stderr.write("analyzed %d mutations in %s!\n" % (len(Ms) ,population))

    # run mutiple regression for bacillus
    # test for interaction and whether there is a difference between treatments
    # create pandas dataframe
    #dict_B_S = {}
    #for population, mutations in mutation_trajectories.items():
    #    #print(population, mutations)
    #    if ('B' in population) or ('S' in population):
    #        dict_B_S[population] = {}
    #        #dict_B_S[population]['taxon'] = population[1]
    #        if population[1] == 'B':
    #            dict_B_S[population]['taxon'] = 0
    #        else:
    #            dict_B_S[population]['taxon'] = 1
    #        dict_B_S[population]['treatment'] = int(population[0])
    #        dict_B_S[population]['replicate'] = int(population[2])
    #        dict_B_S[population]['mutations'] = float(mutations[1])
    #        dict_B_S[population]['genrations_log10'] = np.log10(pt.get_B_S_generations(dict_B_S[population]['taxon'] , str(dict_B_S[population]['treatment']) ))


    #df_B_S = pd.DataFrame.from_dict(dict_B_S).T
    #df_B_S = df_B_S.astype(float)

    #sys.stderr.write("analyzed %d mutations in %s!\n" % (len(Ms) ,population))
    #model_1 = smf.ols(formula='mutations ~ treatment + C(taxon)', data=df_B_S)
    #+ C(taxon)
    #results_1 = model_1.fit()

    #name = ['Lagrange multiplier statistic', 'p-value',
    #    'f-value', 'f p-value']
    #test = sms.het_breuschpagan(results_1.resid, results_1.model.exog)
    #print(lzip(name, test))

    #print(results_1.summary())
    #print(results_1.params)

    #sys.stdout.write("%s\n" % str(res_1.summary()))

    # confidence interval code

    #y_log10_pred = np.asarray([intercept + (slope*x_log_10_i) for x_log_10_i in x_log10])
    #SSE = sum((y_log10 - y_log10_pred) ** 2)
    #N = len(x)
    #sd_SSE = np.sqrt( (1/ (N-2)) * SSE)
    #sxd=np.sum((x_log10-np.mean(x_log10))**2)

    #sx=(x_log10_range-np.mean(x_log10))**2	# x axisr for band
    # Quantile of Student's t distribution for p=1-alpha/2
    #alpha=1-conf
    #q = stats.t.ppf(1-alpha/2, N-2)
    # Confidence band
    #dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
    #print(dy)
    # Upper confidence band
    #ucb = y_log10_range_pred + dy
    # Lower confidence band
    #lcb = y_log10_range_pred - dy

    #ax.plot(10**x_log10_range, 10**lcb, color='k', linestyle=':', linewidth=2)
    #ax.plot(10**x_log10_range, 10**ucb, color='k', linestyle=':', linewidth=2)



    fig = plt.figure(figsize = (14, 12)) #
    fig.subplots_adjust(bottom= 0.15)

    #time_measures = ['transfers', 'generations']

    x_axis_labels = ['Transfer time (days)', 'Generations']
    #slope_null = [-1,1]
    sub_plot_counts = 0

    for taxon_idx, taxon in enumerate(taxa):
        #for time_measure_idx, time_measure in enumerate(time_measures):
        ax_plot = plt.subplot2grid((2, 2), (0, taxon_idx), colspan=1)

        ax_residuals = plt.subplot2grid((2, 2), (1, taxon_idx), colspan=1)

        #if time_measure_idx == 0:
        ax_plot.set_title(latex_dict[taxon], fontsize=24 )

        ax_plot.set_xlabel(x_axis_labels[1], fontsize=20)

        #ax_taxon = plt.subplot2grid((2, 2), tuples[taxon_idx], colspan=1)
        #ax_taxon = fig.add_subplot(2, 2, taxon_idx+1)
        ax_plot.set_xscale('log', basex=10)
        ax_plot.set_yscale('log', basey=10)
        ax_plot.xaxis.set_tick_params(labelsize=16)
        ax_plot.yaxis.set_tick_params(labelsize=16)

        ax_plot.set_ylabel('Mutations at day ' + str(set_time) + ', ' + r'$M({{{}}})$'.format(set_time), fontsize=20)

        ax_plot.text(-0.1, 1.07, sub_plot_labels[sub_plot_counts], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_plot.transAxes)

        #fig.text(0.07, 0.5, 'Mutations at day ' + str(set_time) + ', ' + r'$M({{{}}})$'.format(set_time), ha='center', va='center', rotation='vertical',  fontsize=28)

        ax_residuals.set_xscale('log', basex=10)
        ax_residuals.set_yscale('log', basey=10)
        ax_residuals.xaxis.set_tick_params(labelsize=14)
        ax_residuals.yaxis.set_tick_params(labelsize=14)

        ax_residuals.set_xlabel('Predicted mutations at day ' + str(set_time) + ', ' + r'$\widehat{{M}}({{{}}})$'.format(set_time) , fontsize=18)
        ax_residuals.set_ylabel('Residuals', fontsize=24)
        ax_residuals.text(-0.1, 1.07, sub_plot_labels[sub_plot_counts+2], fontsize=24, fontweight='bold', ha='center', va='center', transform=ax_residuals.transAxes)

        sub_plot_counts += 1

        times_all_list = []
        mutations_all_list = []

        for treatment in treatments:

            mutations_list = np.asarray([value[1] for key, value in mutation_trajectories.items() if treatment+taxon in key])

            if time_measure == 'transfers':
                times_list = np.repeat( int(treatment), len(mutations_list))
            else:
                times_list = np.repeat( np.log10(pt.get_B_S_generations(taxon , treatment)), len(mutations_list))

            ax_plot.scatter((10**times_list) + np.random.randn(len(times_list))*0.1 , 10**mutations_list, s= 180, linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=0.8, zorder=3)
            times_all_list.extend(times_list)
            mutations_all_list.extend(mutations_list)

        ax_plot.set_ylim([(10** min(mutations_all_list))*0.5, (10** max(mutations_all_list))*2  ])
        ax_plot.set_xlim([(10** min(times_all_list))*0.5, (10** max(times_all_list))*2  ])

        slope, intercept, r_value, p_value, std_err = stats.linregress(times_all_list, mutations_all_list)
        x_log10_fit_range =  np.linspace(min(times_all_list) * 0.5, max(times_all_list) * 1.5, 10000)

        y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)
        ax_plot.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', zorder=2)

        y_null_range = 10 ** (slope_null * x_log10_fit_range + intercept)
        ax_plot.plot(10**x_log10_fit_range, y_null_range, c='darkgrey', linestyle=':', lw=3, zorder=1)

        # hypothetical slope of -1
        ratio = (slope - slope_null) / std_err
        pval = stats.t.sf(np.abs(ratio), len(mutations_all_list)-2)*2
        # two sided or one sided?
        sys.stderr.write("Species %s slope t-test = %f, p = %f\n" % (taxon , round(ratio, 3), round(pval, 3)))


        if time_measure == 'transfers':

            ax_plot.text(0.05, 0.22, r'$\beta_{1}=$' + str(round(slope, 2)), fontsize=15, transform=ax_plot.transAxes)
            ax_plot.text(0.05, 0.13, r'$r^{2}=$' + str(round(r_value**2, 2)), fontsize=15, transform=ax_plot.transAxes)

            if pval < 0.05:
                ax_plot.text(0.05, 0.04, r'$\mathrm{p} < 0.05$', fontsize=15, transform=ax_plot.transAxes)
            else:
                ax_plot.text(0.05, 0.04, r'$\mathrm{p} \nless 0.05$', fontsize=15, transform=ax_plot.transAxes)

        else:

            ax_plot.text(0.75, 0.22, r'$\beta_{1}=$' + str(round(slope, 2)), fontsize=15, transform=ax_plot.transAxes)
            ax_plot.text(0.75, 0.13, r'$r^{2}=$' + str(round(r_value**2, 2)), fontsize=15, transform=ax_plot.transAxes)

            if pval < 0.05:
                ax_plot.text(0.75, 0.04, r'$\mathrm{p} < 0.05$', fontsize=15, transform=ax_plot.transAxes)
            else:
                ax_plot.text(0.75, 0.04, r'$\mathrm{p} \nless 0.05$', fontsize=15, transform=ax_plot.transAxes)


        # plot residuals
        y_pred_list = []
        y_resid_list = []

        y_resid_mean = []
        y_pred_mean_list = []
        for treatment in treatments:

            mutations_list = np.asarray([value[1] for key, value in mutation_trajectories.items() if treatment+taxon in key])

            if time_measure == 'transfers':
                times_list = np.repeat( int(treatment), len(mutations_list))
            else:
                times_list = np.repeat( np.log10(pt.get_B_S_generations(taxon , treatment)), len(mutations_list))

            y_pred = (intercept + ( times_list * slope ))

            y_resid = mutations_list - y_pred

            y_pred_list.extend(y_pred)
            y_resid_list.extend(y_resid)

            y_resid_mean.append(np.mean(y_resid))
            y_pred_mean_list.append( intercept + ( np.log10(pt.get_B_S_generations(taxon , treatment)) * slope )  )

            #ax_residuals.scatter((10**y_pred) + np.random.randn(len(times_list))*0.1 , 10**y_resid, s= 110, facecolors=pt.get_colors(treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=0.8, zorder=3)

            #facecolors='none', edgecolors=pt.get_colors(treatment)
            ax_residuals.scatter(10**y_pred, 10**y_resid, s= 180, linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=0.8, zorder=2)

        ax_residuals.set_xlim([(10** min(y_pred_list))*0.5, (10** max(y_pred_list))*2  ])
        #ax_residuals.set_ylim([(10** min(y_resid_list))*0.5, (10** max(y_resid_list))*2  ])

        ax_residuals.set_ylim([ 0.05,  40] )

        ax_residuals.plot(10**np.asarray(y_pred_mean_list), 10**np.asarray(y_resid_mean), c= 'k', marker=pt.plot_species_marker(taxon), linestyle='dashed', ms = 15, lw=2.5, zorder=3)

        ax_residuals.axhline(y=1, color='darkgrey', linestyle=':', lw = 3, zorder=1)



    #fig.text(0.4, 0.04, "Transfer time (days)", va='center', fontsize=28)
    #fig.text(0.07, 0.5, 'Mutations at day ' + str(set_time) + ', ' + r'$M({{{}}})$'.format(set_time), ha='center', va='center', rotation='vertical',  fontsize=28)

    #ax1.text(-0.1, 1.07, "a)", fontsize=11, fontweight='bold', ha='center', va='center', transform=ax1.transAxes)
    #ax2.text(-0.1, 1.07, "b)", fontsize=11, fontweight='bold', ha='center', va='center', transform=ax2.transAxes)
    #ax3.text(-0.1, 1.07, "c)", fontsize=11, fontweight='bold', ha='center', va='center', transform=ax3.transAxes)
    #ax4.text(-0.1, 1.07, "d)", fontsize=11, fontweight='bold', ha='center', va='center', transform=ax4.transAxes)


    fig.subplots_adjust(wspace=0.3) #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/mutation_regression_B_S.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()




def plot_fmax(taxon='S'):

    fig = plt.figure(figsize = (12, 12))

    gene_data = parse_file.parse_gene_list('B')

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    # to get the common gene names for each ID

    alpha_treatment_dict = {'0':0.5, '1':0.5, '2':0.8}


    #ax_mult_freq.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    #ax_mult_freq.set_ylabel('Mean maximum allele frequency, ' + r'$f_{max}$', fontsize=14)

    fmax_treatment_dict = {}

    for treatment_idx, treatment in enumerate(treatments):

        fmax_treatment_dict[treatment] = []

        fmax_list = []

        significan_multiplicity_taxon = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon), "r")
        significan_multiplicity_list = []
        for i, line in enumerate(significan_multiplicity_taxon):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            significan_multiplicity_list.append(items[0])


        populations = [treatment+taxon + replicate for replicate in replicates ]

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))
        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        sys.stdout.write("Total parallelism for %s = %g (p=%g)\n" % (treatment+taxon, G,pvalue))

        for gene_name in convergence_matrix.keys():

            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

            #Ls.append(convergence_matrix[gene_name]['length'])
            #m = gene_parallelism_statistics[gene_name]['multiplicity']

            for population in populations:
                for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:
                    fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                    #predictors.append(m)
                    #responses.append(fixed_weight)

                    #n+=1
                    #nfixed+=fixed_weight

                    #fmax_treatment_dict[treatment].append(f_max)
                    fmax_list.append(f_max)

            #if n > 0.5:
            #    gene_hits.append(n)
            #    gene_predictors.append(m)
            #    #mean_gene_freqs.append(np.mean(freqs))

            #    if nf_max > 0:
            #        ax_mult_freqs_x.append(m)
            #        ax_mult_freqs_y.append( nf_max / n )

        #Ls = np.asarray(Ls)
        #ntot = len(predictors)
        #mavg = ntot*1.0/len(Ls)

        #print(np.mean(fmax_list), np.var(fmax_list))

        plt.hist(fmax_list, label=latex_dict[taxon], linestyle=pt.get_taxon_ls(taxon), color= pt.get_colors(treatment), lw=3, histtype='step', bins = 40, alpha =1, weights=np.zeros_like(fmax_list) + 1. / len(fmax_list))


        # step function
        #ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3.5, color=pt.get_colors(treatment),alpha=0.8, ls=pt.get_taxon_ls(taxon), label= str(int(10**int(treatment))) + '-day', drawstyle='steps', zorder=2)

        #ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color=pt.get_colors(treatment), alpha=0.8, ls=':',  zorder=1)

        #plt.hist(ax_mult_freqs_x, ax_mult_freqs_y, color=pt.get_colors(treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=alpha_treatment_dict[treatment])

        #all_mults.extend(ax_mult_freqs_x)
        #all_freqs.extend(ax_mult_freqs_y)

        #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(ax_mult_freqs_x), np.log10(ax_mult_freqs_y))
        #print(slope, p_value)



    #venn = venn3(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1'), pt.get_colors('2')))
    #c = venn3_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')

    #plt.xscale('log', basex=10)

    fig.suptitle(latex_dict[taxon], fontsize=30)

    fig.subplots_adjust() #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/fmax_%s.png" % taxon
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def likelihood_plot(ntot_subsample=50, subsamples=10000):
    # ntot_subsample minimum number of mutations

    G_subsample_dict = {}

    G_all_mutations_dict = {}

    for taxon in ['B', 'S']:

        for treatment in treatments:

            # Load convergence matrix
            convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

            populations = [treatment+taxon + replicate for replicate in replicates ]

            gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

            G_subsample_list = []
            for i in range(subsamples):

                G_subsample = mutation_spectrum_utils.calculate_subsampled_total_parallelism(gene_parallelism_statistics, ntot_subsample=ntot_subsample)

                G_subsample_list.append(G_subsample)

            G_subsample_list.sort()

            G_subsample_dict[treatment+taxon] = G_subsample_list


            G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

            G_all_mutations_dict[treatment+taxon] = G


    fig = plt.figure(figsize = (8, 8))

    taxon_xaxis_dict = {'B': -0.1, 'S':0.1}

    for treatment in treatments:

        for taxon in ['B', 'S']:

            G_subsample_mean = np.mean(G_subsample_dict[treatment+taxon])
            G_subsample_025 = G_subsample_dict[treatment+taxon][ int( 0.025 * subsamples)  ]
            G_subsample_975 = G_subsample_dict[treatment+taxon][ int( 0.975 * subsamples)  ]

            #xerr1 = [ [z_lclb_mpd_null_mean - lclb_mpd_025, z_lcpl_mpd_null_mean - lcpl_mpd_025, z_hclb_mpd_null_mean - hclb_mpd_025, z_hcpl_mpd_null_mean - hcpl_mpd_025 ] ,
            #        [lclb_mpd_975 - z_lclb_mpd_null_mean, lcpl_mpd_975 - z_lcpl_mpd_null_mean, hclb_mpd_975 - z_hclb_mpd_null_mean, hcpl_mpd_975 -z_hcpl_mpd_null_mean ]]

            plt.errorbar(int(treatment) + taxon_xaxis_dict[taxon], G_subsample_mean, yerr = [ [G_subsample_mean-G_subsample_025], [ G_subsample_975-G_subsample_mean]], \
                    fmt = 'o', alpha = 1, barsabove = True, marker = 's', \
                    mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=1, ms=17)

            plt.scatter(int(treatment) + taxon_xaxis_dict[taxon], G_subsample_mean, marker='s', s = 250, \
                linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)


            plt.scatter(int(treatment) + taxon_xaxis_dict[taxon], G_all_mutations_dict[treatment+taxon], marker=pt.plot_species_marker(taxon), s = 250, \
                linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)

    plt.xlabel("Transfer time (days)", fontsize = 20)

    plt.xticks((0,1,2), ('1', '10', '100'), fontsize=14  )
    plt.rc('ytick',labelsize=12)

    plt.ylim([1.2,6.2])


    plt.ylabel("Net increase in log-likelihood, " r'$\Delta \ell$' , fontsize = 20)



    legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=latex_dict['B'],
                        markerfacecolor='k', markersize=13),
                    Line2D([0], [0], marker='o', color='none', label=latex_dict['S'],
                        markerfacecolor='w', markersize=13, markeredgewidth=2),
                    Line2D([0], [0], color = 'none', marker='s', label=latex_dict['B'] + ', sub-sampled',
                        markerfacecolor='k', markersize=13),
                    Line2D([0], [0], marker='s', color='none', label=latex_dict['S'] + ', sub-sampled',
                        markerfacecolor='w', markersize=13, markeredgewidth=2)]
    # Create the figure
    plt.legend(handles=legend_elements, loc='upper left')

    fig.subplots_adjust() #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/G_score_subsample.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    Returns
    -------
    matplotlib.patches.Ellipse
    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)



def plot_pca_B_S():

    pop_mutations_dict = {}

    for treatment in treatments:

        for taxon in ['B', 'S']:

            populations = [treatment+taxon + replicate for replicate in replicates ]

            # Load convergence matrix
            convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))
            gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations, Lmin=100)
            mean_length = np.mean([gene_parallelism_statistics[key]['length'] for key in gene_parallelism_statistics.keys() ])
            for gene_name in convergence_matrix.keys():

                # get multiplicities for regression
                if gene_name not in pop_mutations_dict:
                    pop_mutations_dict[gene_name] = {}

                convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

                for population in populations:

                    fixed_weight_sum = 0

                    for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:

                        fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                        fixed_weight_sum += fixed_weight

                    pop_mutations_dict[gene_name][population] = fixed_weight_sum * (mean_length / gene_parallelism_statistics[gene_name]['length'] )



    df = pd.DataFrame.from_dict(pop_mutations_dict)
    df = df.fillna(0)
    df = df.loc[:, (df != 0).any(axis=0)]
    df_np = df.values
    X = df_np/df_np.sum(axis=1)[:,None]
    #X -= np.mean(X, axis = 0)

    pca = PCA()
    df_out = pca.fit_transform(X)

    fig = plt.figure(figsize = (6, 8))
    fig.tight_layout(pad = 2.8)
    # Scatterplot on main ax
    ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1)
    ax1.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
    ax1.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=2)
    ax1.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=3)
    ax1.set_xlabel('PC 1 (' + str(round(pca.explained_variance_ratio_[0]*100,2)) + '%)' , fontsize = 12)
    ax1.set_ylabel('PC 2 (' + str(round(pca.explained_variance_ratio_[1]*100,2)) + '%)' , fontsize = 12)


    df_out_labelled = pd.DataFrame(data=df_out, index=df.index.values)

    for treatment in treatments:

        for taxon in ['B', 'S']:

            #for replicate in replicates:
            treatment_taxon_idx_list = [i for i  in df_out_labelled.index if treatment+taxon in i]
            df_treatment_taxon = df_out_labelled.loc[ treatment_taxon_idx_list , : ]

            ax1.scatter(df_treatment_taxon.values[:,0], df_treatment_taxon.values[:,1], marker=pt.plot_species_marker(taxon), facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), linewidth=2, alpha = 0.8, s = 220, zorder=4)

            confidence_ellipse(df_treatment_taxon.ix[:,0],df_treatment_taxon.ix[:,1], ax1, n_std=2, edgecolor=pt.get_colors(treatment), linestyle=pt.get_taxon_ls(taxon), lw=3)

    ax1.set_xlim([-1,1])
    ax1.set_ylim([-1,1])


    #ax1.text(0.7,0.8,r'$n_{\mathrm{ESCRE1901}}=$' + str(len(df_pops_ESCRE1901)), fontsize=11, color='r', ha='center', va='center', transform=ax1.transAxes  )
    #ax1.text(0.7,0.7,r'$n_{\mathrm{ECB \_ 01992 }} = $' + str(len(df_pops_ECB_01992)), fontsize=11, color='purple', ha='center', va='center', transform=ax1.transAxes)


    plt.tight_layout()
    fig.savefig(pt.get_path() + '/figs/pca.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






#likelihood_plot()


#plot_bPTR_all()

#allele_survival()


#get_mult_similarity()

#plot_fmax()

#plot_spo0a_fitness()
#plot_spores()

#temporal_plasmid_coverage_B_S()
#temporal_plasmid_coverage_D()
#plot_mutation_trajectory_B_S()
#plot_allele_corr_delta_B_S()
#mutation_regression_B_S()

#taxa = ['B', 'S', 'C', 'D', 'F', 'J', 'P']
#for taxon in taxa:
#    #plot_within_taxon_paralleliism(taxon)
#    if (taxon != 'B') and (taxon != 'S'):
#        plot_mutation_trajectory_taxon(taxon)



#plot_B_S_paralleliism()


#plot_allele_corr_delta()
#plot_allele_corr_delta_B_S()
