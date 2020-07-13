from __future__ import division
import os, sys
import numpy as np
import mutation_spectrum_utils
import phylo_tools as pt

import  matplotlib.pyplot as plt
import statsmodels.stats.multitest as mt


treatments=pt.treatments
replicates = pt.replicates
iter = 100000

taxa = ['B', 'C', 'D', 'F', 'J', 'P']

extinction_counts = {}

file = open(pt.get_path() + '/data/Tranfers_history2016.txt')
for line in file:
    line = line.strip().split('\t')
    population = line[0]
    cause = line[5]
    if 'did not grow' in cause:
        population_split = population.split('-')

        population_id = population_split[0] + population_split[-2] + population_split[-1]

        if population_id not in extinction_counts:
            extinction_counts[population_id] = 0

        extinction_counts[population_id] += 1

file.close()

CV_params = []

for taxon in taxa:
    for treatment in treatments:

        taxon_treatment_pops = [extinction_counts[i] for i in list(extinction_counts.keys()) if treatment+taxon in i]

        if len(taxon_treatment_pops) < 3:
            continue

        full_counts = np.asarray(taxon_treatment_pops + [0] * (5-len(taxon_treatment_pops)))

        mean_counts = np.mean(full_counts)

        cv_obs = np.std(full_counts) / mean_counts

        cv_null = []
        for j in range(iter):
            counts_pois = np.random.poisson(lam=mean_counts, size=len(full_counts))
            cv_null_j = np.std(counts_pois) / np.mean(counts_pois)
            # lambda ** (-0.5) = 0 if lambda = 0
            if np.isnan(cv_null_j) == True:
                cv_null_j = 0
            cv_null.append( cv_null_j)

        cv_null.sort()
        cv_null = np.asarray(cv_null)

        p_value = len([x for x in cv_null if x > cv_obs]) / iter

        #cv_null_standard = (cv_null - np.mean(cv_null)) / np.std(cv_null)
        #cv_obs_standard = (cv_obs - np.mean(cv_null)) / np.std(cv_null)
        #CI_025 = cv_null_standard[int(10000*0.025)]
        #CI_975 = cv_null_standard[int(10000*0.975)]

        CV_params.append([treatment+taxon, mean_counts, cv_obs, p_value])


p_values = np.asarray([CV_param[-1] for CV_param in CV_params])

sig_p_values_idx = np.where(p_values < 0.5)[0]

reject, pvals_corrected, alphacSidak, alphacBonf = mt.multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)


for CV_param_idx, CV_params in enumerate(CV_params):
    sys.stdout.write("%s n_ext=%f, CV=%f, p=%f \n" % (CV_params[0], CV_params[1], CV_params[2], pvals_corrected[CV_param_idx] ) )



#fig, ax = plt.subplots(figsize=(8,8))

#x_axis_tick_labels = [""]

#for CV_param_idx, CV_param in enumerate(CV_params):
#    taxon = CV_param[0][1]
#    treatment = CV_param[0][0]
#    print(taxon) # print op value and CV

#    x_axis_tick_labels.append(pt.get_treatment_name(treatment))

#    ax.errorbar(CV_param_idx, 0, yerr = [ [0-CV_param[2]], [ CV_param[3]-0]], \
#            fmt = 'o', alpha = 1, barsabove = True, marker = '.', \
#            mfc = 'k', mec = 'k', c = 'k', lw=3.5, ms=30, zorder=2)

#    #plt.scatter(CV_param_idx, 0, marker='s', s = 250, \
#    #    linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)

#    ax.scatter(CV_param_idx, CV_param[1], marker=pt.plot_species_marker(taxon), s = 200, \
#        linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), \
#        edgecolors=pt.get_colors(treatment), alpha=1, zorder=3)


#ax.axhline(y=0, color='grey', linestyle=':', lw=2,alpha = 0.8, zorder=1)
#ax.set_ylabel("Standardized CV of extinction events", fontsize = 20)
#ax.set_xticklabels(x_axis_tick_labels)

#ax.text(0.15, -0.1, r'$\mathbf{\mathit{Caulobacter}}$', fontsize=16, fontweight='bold', ha='center', va='center', transform=ax.transAxes)
#ax.text(0.5, -0.1, r'$\mathbf{\mathit{Deinococcus}}$', fontsize=16, fontweight='bold', ha='center', va='center', transform=ax.transAxes)
#ax.text(0.85, -0.1, r'$\mathbf{\mathit{Janthinobacterium}}$', fontsize=16, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

#fig.subplots_adjust() #hspace=0.3, wspace=0.5
#fig_name = pt.get_path() + "/figs/extinctions.pdf"
#fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.close()
