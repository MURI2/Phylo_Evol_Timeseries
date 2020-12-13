from __future__ import division
import os, sys, itertools
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy as np
import pandas as pd
from itertools import combinations

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import scipy.stats as stats
import statsmodels.stats.multitest as multitest

np.random.seed(123456789)


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i', 'j','k','l']
all_subplot_counts = 0

n_permutations = 10000
#treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']

#taxa = ['B']

legend_elements = [Line2D([0], [0], color=pt.get_colors('0'), lw=2, label='1-day'),
                   Line2D([0], [0], color=pt.get_colors('1'), lw=2, label='10-days'),
                  Line2D([0], [0], color=pt.get_colors('2'), lw=2, label='100-day')]


nonsynonymous_types = set(['missense','nonsense'])
synonymous_types = set(['synonymous'])

non_appeared = {}
non_fixed = {}

syn_appeared = {}
syn_fixed = {}

targeted_Lsyn = {}
targeted_Lnon = {}
targeted_fixed_Lsyn = {}
targeted_fixed_Lnon = {}

taxon_Lsyn_dict = {}
taxon_Lnon_dict = {}

populations = []

fmax_dict = {}

for taxon in taxa:
    Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes(taxon)
    taxon_Lsyn_dict[taxon] = Lsyn
    taxon_Lnon_dict[taxon] = Lnon
    fmax_dict[taxon] = {}

    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment in treatments:

        nonsynonymous_fmax_all = []
        synonymous_fmax_all = []

        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            populations.append(population)

            non_appeared[population] = 1
            non_fixed[population] = 1

            syn_appeared[population] = 1
            syn_fixed[population] = 1

            targeted_Lsyn[population] = 1
            targeted_Lnon[population] = 1
            targeted_fixed_Lsyn[population] = 1
            targeted_fixed_Lnon[population] = 1

            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            num_processed_mutations = 0

            for mutation_idx in range(0,len(mutations)):

                #location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                location, gene_name, allele, var_type, codon, position_in_codon, AAs_count,  test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

                state_Ls = state_trajectories[mutation_idx]

                good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]
                masked_state_Ls = state_Ls[good_idxs]

                fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1],masked_freqs[-1])

                if var_type in nonsynonymous_types or var_type in synonymous_types:
                    targeted_Lnon[population] += (1-substitution_specific_synonymous_fraction[allele])
                    targeted_fixed_Lnon[population] += fixed_weight*(1-substitution_specific_synonymous_fraction[allele])
                    targeted_Lsyn[population] += substitution_specific_synonymous_fraction[allele]
                    targeted_fixed_Lsyn[population] += fixed_weight*substitution_specific_synonymous_fraction[allele]

                if var_type in nonsynonymous_types:
                    non_appeared[population]+=1
                    non_fixed[population]+=fixed_weight
                    num_processed_mutations+=1

                    nonsynonymous_fmax_all.append(max(freqs))

                elif var_type in synonymous_types:
                    syn_appeared[population]+=1
                    syn_fixed[population]+=fixed_weight
                    num_processed_mutations+=1

                    synonymous_fmax_all.append(max(freqs))

        fmax_dict[taxon][treatment] = {}
        fmax_dict[taxon][treatment]['synonymous_fmax'] = np.asarray(synonymous_fmax_all)
        fmax_dict[taxon][treatment]['nonsynonymous_fmax'] = np.asarray(nonsynonymous_fmax_all)



#fig, ax = plt.subplots(figsize=(4,4))
fig = plt.figure(figsize = (9, 6))
gs = gridspec.GridSpec(nrows=2, ncols=3)
dn_ds_count = 0
all_subplot_counts = 0
fmax_min = 0.1
fmax_max = 0.42
fmax_max_range = np.linspace(fmax_min, fmax_max, num=1000)
fmax_cutoff_dict = {}
#for taxon in taxa:
for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):
    for taxon_idx, taxon in enumerate(taxon_list):

        ax = fig.add_subplot(gs[taxon_list_idx, taxon_idx])
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )
        ax.text(-0.1, 1.07, sub_plot_labels[all_subplot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

        if all_subplot_counts == 0:
            ax.legend(handles=legend_elements, loc='lower left')

        all_subplot_counts+=1

        fmax_dict_taxon = fmax_dict[taxon]

        fmax_cutoff_dict[taxon] = {}

        for fmax_dict_taxon_key in fmax_dict_taxon.keys():

            fmax_dict_taxon_treatment_dict = fmax_dict_taxon[fmax_dict_taxon_key]

            nonsynonymous_fmax = fmax_dict_taxon_treatment_dict['nonsynonymous_fmax']
            synonymous_fmax = fmax_dict_taxon_treatment_dict['synonymous_fmax']

            pNpS_fmax = []

            fmax_max_passed = []

            for fmax_max in fmax_max_range:

                nonsynonymous_fmax_max = nonsynonymous_fmax[nonsynonymous_fmax>=fmax_max]
                synonymous_fmax_max = synonymous_fmax[synonymous_fmax>=fmax_max]

                n_nonsynonymous_fmax_max = len(nonsynonymous_fmax_max) + 1
                n_synonymous_fmax_max = len(synonymous_fmax_max) + 1

                if (n_nonsynonymous_fmax_max == 0) or (n_synonymous_fmax_max == 0):
                    continue

                pNpS = (n_nonsynonymous_fmax_max/n_synonymous_fmax_max) * ((taxon_Lsyn_dict[taxon]+1)/ (taxon_Lnon_dict[taxon]+1))

                fmax_max_passed.append(fmax_max)
                pNpS_fmax.append(pNpS)

            pNpS_fmax = np.asarray(pNpS_fmax)
            fmax_cutoff_dict[taxon][fmax_dict_taxon_key] = pNpS_fmax

            ax.plot(fmax_max_passed, pNpS_fmax, c =pt.get_colors(fmax_dict_taxon_key), linewidth=2, zorder=2)

            ax.set_xlim([fmax_min, fmax_max])

            if len(pNpS_fmax[pNpS_fmax>=1]) != 0:
                ax.axhline(y=1, color='k', linestyle=':', lw=2.5, alpha = 0.8, zorder=1)



fig.text(0.5, -0.01, 'Maximum observed allele frequency, ' + r'$f_{max}$', ha='center', va='center', fontsize=18)
fig.text(-0.01, 0.5, r'$pN/pS$' + ' for mutations ' + r'$\geq f_{max}$', ha='center', va='center', rotation='vertical', fontsize=18)

fig.subplots_adjust(hspace=0.4, wspace=0.6) #hspace=0.3, wspace=0.5
fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/dn_ds_fmax.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()





record_strs = [",".join(['treatment_pair', 'taxon', 'tree_name', 'mean_absolute_difference'])]

msd_dict = {}
for taxon in taxa:

    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment_pair in combinations(treatments, 2):

        if treatment_pair not in msd_dict:
            msd_dict[treatment_pair] = {}

        sample_1 = fmax_cutoff_dict[taxon][treatment_pair[0]]
        sample_2 = fmax_cutoff_dict[taxon][treatment_pair[1]]

        mean_absolute_difference = np.mean(np.absolute(sample_1 - sample_2))

        msd_dict[treatment_pair][taxon] = {}
        msd_dict[treatment_pair][taxon]['mean_absolute_difference'] = mean_absolute_difference


        record_str = ",".join(['%s_%s' % treatment_pair,  str(taxon), pt.tree_name_dict[taxon], str(mean_absolute_difference)])
        record_strs.append(record_str)




sys.stderr.write("Writing intermediate file for phylogenetic ANOVA...\n")
file = open(pt.get_path()+'/data/pNpS_mean_absolute_differences.csv',"w")
record_str = "\n".join(record_strs)
file.write(record_str)
file.close()
sys.stderr.write("Done!\n")


for treatment_pair in combinations(pt.treatments, 2):

    mean_absolute_difference_list = [msd_dict[treatment_pair][taxon_i]['mean_absolute_difference'] for taxon_i in msd_dict[treatment_pair].keys()]

    mean_absolute_difference_mean = np.mean(mean_absolute_difference_list)
    mean_absolute_difference_se =  np.std(mean_absolute_difference_list) / np.sqrt(len(mean_absolute_difference_list))

    sys.stderr.write("%d vs. %d Mean MAD = %f, SE MAD = %f\n" % (10**int(treatment_pair[0]), 10**int(treatment_pair[1]), mean_absolute_difference_mean, mean_absolute_difference_se))





taxa_to_test = ['B','C','D','F','P']

within_sum = 0
all_divergenes_to_test = []
all_mean_divergenes = []
all_n_divergences = []
for treatment_pair in msd_dict.keys():

    divergences_treatment_pair = [msd_dict[treatment_pair][taxon]['mean_absolute_difference'] for taxon in taxa_to_test]
    all_divergenes_to_test.extend(divergences_treatment_pair)

    divergences_treatment_pair = np.asarray(divergences_treatment_pair)

    within_sum += sum((divergences_treatment_pair - np.mean(divergences_treatment_pair))**2)

    all_mean_divergenes.append(np.mean(divergences_treatment_pair))
    all_n_divergences.append(len(divergences_treatment_pair))


all_divergenes_to_test = np.asarray(all_divergenes_to_test)
all_mean_divergenes = np.asarray(all_mean_divergenes)
all_n_divergences = np.asarray(all_n_divergences)
F_numerator = sum(all_n_divergences*((all_mean_divergenes - np.mean(all_divergenes_to_test))**2)) /( len(all_mean_divergenes)-1 )

F_denominator = within_sum / (len(all_divergenes_to_test) - len(all_mean_divergenes))
F = F_numerator/F_denominator


F_permute_all = []
for i in range(n_permutations):

    vs_1_10 = []
    vs_1_100 = []
    vs_10_100 = []

    for taxon in taxa_to_test:

        taxon_standardized_corr = []
        treatment_pairs_list = msd_dict.keys()
        divergences_taxon = np.asarray([msd_dict[l][taxon]['mean_absolute_difference'] for l in treatment_pairs_list])
        divergences_taxon_permute = np.random.permutation(divergences_taxon)
        vs_1_10.append(divergences_taxon_permute[0])
        vs_1_100.append(divergences_taxon_permute[1])
        vs_10_100.append(divergences_taxon_permute[2])

    vs_1_10 = np.asarray(vs_1_10)
    vs_1_100 = np.asarray(vs_1_100)
    vs_10_100 = np.asarray(vs_10_100)

    vs_all = np.concatenate((vs_1_10, vs_1_100, vs_10_100))
    vs_arrays = [vs_1_10, vs_1_100, vs_10_100]
    within_sum_permute = 0
    all_means_permute = []
    all_n_divergences_permute = []
    for vs_array in vs_arrays:
        within_sum_permute += sum((vs_array - np.mean(vs_array))**2)
        all_means_permute.append(np.mean(vs_array))
        all_n_divergences_permute.append(len(vs_array))

    means_all_permute = np.asarray(all_means_permute)

    F_permute_numerator = sum((all_n_divergences_permute*((all_means_permute - np.mean(all_means_permute))**2)))/ (len(all_means_permute)-1)
    F_permute_denominator = within_sum_permute/ (len(vs_all) - len(all_means_permute))

    F_permute = F_permute_numerator/F_permute_denominator
    F_permute_all.append(F_permute)

F_permute_all = np.asarray(F_permute_all)
P_F = len(F_permute_all[F_permute_all>F]) / n_permutations
sys.stderr.write("Mean MAD: F = %.3f, P = %.4f \n" % (F, P_F))






def get_ks_distance():

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
                ks_dict[treatment_pair]['D'] = []


            sample_1 = fmax_cutoff_dict[taxon][treatment_pair[0]]
            sample_2 = fmax_cutoff_dict[taxon][treatment_pair[1]]

            D, p_value = stats.ks_2samp(sample_1, sample_2)

            ks_dict[treatment_pair][taxon] = {}
            ks_dict[treatment_pair][taxon]['D'] = D
            #ks_dict[treatment_pair]['D'] .append() = D
            #ks_dict[treatment_pair][taxon]['p_value'] = p_value

            #ks_dict[treatment_pair]['D'].append(D)




    taxa_to_test = ['B','C','D','F','P']

    within_sum = 0
    all_divergenes_to_test = []
    all_mean_divergenes = []
    all_n_divergences = []
    for treatment_pair in ks_dict.keys():

        divergences_treatment_pair = [ks_dict[treatment_pair][taxon]['D'] for taxon in taxa_to_test]
        all_divergenes_to_test.extend(divergences_treatment_pair)

        divergences_treatment_pair = np.asarray(divergences_treatment_pair)

        within_sum += sum((divergences_treatment_pair - np.mean(divergences_treatment_pair))**2)

        all_mean_divergenes.append(np.mean(divergences_treatment_pair))
        all_n_divergences.append(len(divergences_treatment_pair))


    all_divergenes_to_test = np.asarray(all_divergenes_to_test)
    all_mean_divergenes = np.asarray(all_mean_divergenes)
    all_n_divergences = np.asarray(all_n_divergences)
    F_numerator = sum(all_n_divergences*((all_mean_divergenes - np.mean(all_divergenes_to_test))**2)) /( len(all_mean_divergenes)-1 )

    F_denominator = within_sum / (len(all_divergenes_to_test) - len(all_mean_divergenes))
    F = F_numerator/F_denominator
    # write code to permute while controlling for taxon identiy
    # hacky, but gets the job done
    F_permute_all = []
    for i in range(n_permutations):

        vs_1_10 = []
        vs_1_100 = []
        vs_10_100 = []

        for taxon in taxa_to_test:

            taxon_standardized_corr = []
            treatment_pairs_list = ks_dict.keys()
            divergences_taxon = np.asarray([ks_dict[l][taxon]['D'] for l in treatment_pairs_list])
            divergences_taxon_permute = np.random.permutation(divergences_taxon)
            vs_1_10.append(divergences_taxon_permute[0])
            vs_1_100.append(divergences_taxon_permute[1])
            vs_10_100.append(divergences_taxon_permute[2])

        vs_1_10 = np.asarray(vs_1_10)
        vs_1_100 = np.asarray(vs_1_100)
        vs_10_100 = np.asarray(vs_10_100)

        vs_all = np.concatenate((vs_1_10, vs_1_100, vs_10_100))
        vs_arrays = [vs_1_10, vs_1_100, vs_10_100]
        within_sum_permute = 0
        all_means_permute = []
        all_n_divergences_permute = []
        for vs_array in vs_arrays:
            within_sum_permute += sum((vs_array - np.mean(vs_array))**2)
            all_means_permute.append(np.mean(vs_array))
            all_n_divergences_permute.append(len(vs_array))

        means_all_permute = np.asarray(all_means_permute)

        F_permute_numerator = sum((all_n_divergences_permute*((all_means_permute - np.mean(all_means_permute))**2)))/ (len(all_means_permute)-1)
        F_permute_denominator = within_sum_permute/ (len(vs_all) - len(all_means_permute))

        F_permute = F_permute_numerator/F_permute_denominator
        F_permute_all.append(F_permute)

    F_permute_all = np.asarray(F_permute_all)
    P_F = len(F_permute_all[F_permute_all>F]) / n_permutations
    sys.stderr.write("Mean divergence: F = %.3f, P = %.4f \n" % (F, P_F))




#D_mean_dict = {}
#D_all_samples = []
#D_means = []
#D_n = []

#within_difference_squared = 0

#for treatment_pair in combinations(pt.treatments, 2):

#    D_mean = np.mean(ks_dict[treatment_pair]['D'])
#    D_all_samples.extend(ks_dict[treatment_pair]['D'])
#    D_se =  np.std(ks_dict[treatment_pair]['D']) / np.sqrt(len(ks_dict[treatment_pair]['D']))

#    sys.stderr.write("%d vs. %d Mean D = %f, SE D = %f\n" % (10**int(treatment_pair[0]), 10**int(treatment_pair[1]), D_mean, D_se))

    #D_mean_dict[treatment_pair] = D_mean
#    D_means.append(D_mean)
#    D_n.append(len(ks_dict[treatment_pair]['D']))


#    within_difference_squared += sum((ks_dict[treatment_pair]['D'] - D_mean)**2)


#D_n = np.asarray(D_n)
##D_all_samples = np.asarray(D_all_samples)
#D_overall_mean = np.mean(D_all_samples)
#D_means = np.asarray(D_means)


# calculate F statistic
# between-group variability
#between_variability = sum(D_n*((D_means - D_overall_mean)**2) ) / (len(D_means)-1)

# within-group variability
#within_variability = within_difference_squared / (len(D_all_samples) - len(D_means))

#F = between_variability/within_variability

#F_permute_list = []


#n_permutations= 10000

#def split_array(array, sizes):

#    list_of_arrays = []
#    n_ = 0
#    for size in sizes:
#        list_of_arrays.append(array[n_: n_+size])
#        n_ += size
#    return list_of_arrays


#for i in range(n_permutations):

#    D_all_samples_permute = np.random.permutation(D_all_samples)
    #D_all_treatment_permute = np.split(D_all_samples_permute, D_n)
#    D_all_treatment_permute = split_array(D_all_samples_permute, D_n)
#    within_difference_squared_permute = 0
#    between_difference_squared_permute = 0
#    for D_all_treatment_permute_i_idx, D_all_treatment_permute_i in enumerate(D_all_treatment_permute):

#        D_mean_permute_i = np.mean(D_all_treatment_permute_i)

#        within_difference_squared_permute += sum((D_all_treatment_permute_i - D_mean_permute_i)**2)

#        between_difference_squared_permute += D_n[D_all_treatment_permute_i_idx]*((D_mean_permute_i - D_overall_mean)**2)

#    within_variability_permute = within_difference_squared_permute/(len(D_all_samples) - len(D_means))

#    between_variability_permute = between_difference_squared_permute / (len(D_means)-1)

#    F_permute_list.append(within_variability_permute / between_variability_permute)


#F_permute_list = np.asarray(F_permute_list)


#P = len(F_permute_list[F_permute_list>F]) / n_permutations

#sys.stderr.write("Permutational ANOVA F = %f, P = %f\n" % (F, P))
