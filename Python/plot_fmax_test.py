from __future__ import division
import os, sys, itertools
import numpy as np
from itertools import combinations

import phylo_tools as pt
import parse_file
import timecourse_utils

import  matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import scipy.stats as stats
import statsmodels.stats.multitest as multitest


np.random.seed(123456789)

fmax_dict = {}
#taxa = pt.taxa
taxa = ['B','C','D','P', 'F', 'J']
# loop through taxa and get M(700) for all reps in each treatment
for treatment in pt.treatments:

    fmax_dict[treatment] = {}

for taxon in taxa:

    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment in treatments:

        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        f_max_all = []

        #for population in populations:
        for replicate in pt.replicates:
            population = treatment + taxon + replicate

            if population in pt.populations_to_ignore:
                continue

            for gene_name in sorted(convergence_matrix.keys()):

                for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:

                    f_max_all.append(f_max)

        fmax_dict[treatment][taxon] =  np.asarray(f_max_all)


ks_dict = {}
treatment_pairs = []
p_values = []
for taxon in taxa:

    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment_pair in combinations(treatments, 2):

        if treatment_pair not in ks_dict:
            ks_dict[treatment_pair] = {}

        sample_1 = fmax_dict[treatment_pair[0]][taxon]
        sample_2 = fmax_dict[treatment_pair[1]][taxon]

        D, p_value = stats.ks_2samp(sample_1, sample_2)

        ks_dict[treatment_pair][taxon] = {}
        ks_dict[treatment_pair][taxon]['D'] = D
        ks_dict[treatment_pair][taxon]['p_value'] = p_value

        treatment_pairs.append((treatment_pair, taxon))
        p_values.append(p_value)


reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')
for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):
    ks_dict[treatment_pair[0]][treatment_pair[1]]['p_value_bh'] = pvals_corrected[treatment_pair_idx]


fig, ax = plt.subplots(figsize=(4,4))

for treatment in pt.treatments:

    for taxon, f_max_array in fmax_dict[treatment].items():

        f_max_array_sort = np.sort(f_max_array)
        cdf = 1-  np.arange(len(f_max_array_sort))/float(len(f_max_array_sort))


        #num_bins = 40
        #counts, bin_edges = np.histogram(f_max_array_sort, bins=num_bins, normed=True)
        #cdf = np.cumsum(counts)
        #pylab.plot(bin_edges[1:], cdf)

        ax.plot(f_max_array_sort, cdf, c =pt.get_colors(treatment), alpha=0.8)
        #marker=pt.plot_species_marker(taxon), markersize=1)

        #print(np.linspace(0, len()))



ax.set_xlim([ 0.09, 1.03 ])
ax.set_ylim([ 0.0008, 1.03 ])

ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)

ax.set_xlabel('Maximum observed allele frequency, ' + r'$f_{max}$', fontsize=12)
ax.set_ylabel('Fraction of mutations ' + r'$\geq f_{max}$', fontsize=13)


ins_ks = inset_axes(ax, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.18,0.05,0.4,0.38), bbox_transform=ax.transAxes)
#ins_ks.set_xlabel('Max.' + r'$\left \langle x(t) \right \rangle$', fontsize=8)
ins_ks.set_ylabel("Mean KS distance", fontsize=8)


def rand_jitter(arr):
    stdev = 5 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev


count = 0
mean_D_dict = {}
for treatment_pair_idx, treatment_pair in enumerate(ks_dict.keys()):
    x = []
    y = []
    markers = []
    for treatment_pair_taxon in ks_dict[treatment_pair]:

        D_i = ks_dict[treatment_pair][treatment_pair_taxon]['D']
        p_value_bh_i = ks_dict[treatment_pair][treatment_pair_taxon]['p_value_bh']

        x.append(count)
        y.append(D_i)
        markers.append(pt.plot_species_marker(treatment_pair_taxon))

        if p_value_bh_i < 0.05:
            facecolor_i = 'none'
        else:
            facecolor_i = 'k'

        count_jitter = count +( np.random.randn() * 0.15)

        #marker_style = dict(color='k', marker=pt.plot_species_marker(treatment_pair_taxon),
        #                    markerfacecoloralt=pt.get_colors(treatment_pair[0]),
        #                    markerfacecolor=pt.get_colors(treatment_pair[1]))


        #ins_ks.plot(count_jitter, D_i, markersize = 6.5,   \
        #    linewidth=0.2,  alpha=0.8, fillstyle='left', **marker_style)


    mean_D_dict[treatment_pair] = {}
    mean_D = np.mean(y)
    se_D = np.std(y) / np.sqrt(len(y))
    mean_D_dict[treatment_pair]['mean_D'] = mean_D
    mean_D_dict[treatment_pair]['se_D'] = se_D

    marker_style = dict(color='k',
                        markerfacecoloralt=pt.get_colors(treatment_pair[0]),
                        markerfacecolor=pt.get_colors(treatment_pair[1]))

    if treatment_pair_idx == 0:
        count_real = count + 0.25
    elif treatment_pair_idx == 1:
        count_real = count
    else:
        count_real = count - 0.25

    ins_ks.errorbar(count_real, mean_D, yerr = se_D , \
            fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
            mfc = 'white', mec = 'white', lw=1.5, c = 'k', zorder=1, ms=17)

    ins_ks.plot(count_real, mean_D, markersize = 8, marker = 'o',  \
        linewidth=0.2,  alpha=1, fillstyle='left', zorder=2, **marker_style)


    count+=1







ins_ks.tick_params(labelsize=5)
ins_ks.tick_params(axis='both', which='major', pad=1)

ins_ks.set_xticks([0.15, 1, 1.85])
ins_ks.set_xticklabels(['1 vs. 10-days', '1 vs. 100-days', '10 vs. 100-days'],fontweight='bold' )
ins_ks.tick_params(axis='x', labelsize=3.3, length = 0)

temp = ins_ks.yaxis.get_ticklabels()
temp = list(set(temp) - set(temp[::2]))
for label in temp:
    label.set_visible(False)

from matplotlib.ticker import FormatStrFormatter

ins_ks.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


ins_ks.set_xlim([-0.2,2.2])
ins_ks.set_ylim([0.05,0.25])

#ins_ks.axvline(x=0.5, color='k', linestyle='-', lw=0.9, alpha = 1, zorder=1)
#ins_ks.axvline(x=1.5, color='k', linestyle='-', lw=0.9, alpha = 1, zorder=1)

#ins_ks.axhline(y=0, color='k', linestyle=':', alpha = 1, zorder=1)


fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/fmax_text.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
