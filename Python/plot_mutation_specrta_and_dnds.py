from __future__ import division
import os, sys, itertools
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy as np
import pandas as pd

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

iter=10
iter_corr= 10

spectra_latex_dict = {'AT_GC':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{G}:\mathrm{C}$',
                    'AT_CG':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{C}:\mathrm{G}$',
                    'AT_TA':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_AT':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{A}:\mathrm{T}$',
                    'GC_TA':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_CG':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{C}:\mathrm{G}$'}


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i', 'j','k','l']
all_subplot_counts = 0

treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']
iter
spectra_dict = {}
snps_indel_dict = {}
for taxon in taxa:

    for treatment in treatments:

        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            mutation_counts = {'GC_AT': 1,
                            'GC_TA': 1,
                            'GC_CG': 1,
                            'AT_GC': 1,
                            'AT_CG': 1,
                            'AT_TA': 1}

            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            for mutation_idx in range(0,len(mutations)):

                location, gene_name, allele, var_type, codon, position_in_codon, AAs_count,  test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                state_Ls = state_trajectories[mutation_idx]

                if AAs_count == 1:

                    ref, alt = allele.split('->')

                    if (ref=='A' and alt=='T') or (ref=='T' and alt=='A'):
                        mutation_counts['AT_TA']+=1

                    elif (ref=='A' and alt=='G') or (ref=='T' and alt=='C'):
                        mutation_counts['AT_GC']+=1

                    elif (ref=='A' and alt=='C') or (ref=='T' and alt=='G'):
                        mutation_counts['AT_CG']+=1

                    elif (ref=='G' and alt=='A') or (ref=='C' and alt=='T'):
                        mutation_counts['GC_AT']+=1

                    elif (ref=='G' and alt=='T') or (ref=='C' and alt=='A'):
                        mutation_counts['GC_TA']+=1

                    elif (ref=='G' and alt=='C') or (ref=='C' and alt=='G'):
                        mutation_counts['GC_CG']+=1

                    else:

                        print(ref, alt)

            mutation_counts_rel = {k: v / sum(mutation_counts.values()) for k, v in mutation_counts.items()}

            spectra_dict[population] = mutation_counts_rel


# get indel data now






df = pd.DataFrame.from_dict(spectra_dict)
X = df.values.transpose()

X_std=StandardScaler().fit_transform(X) # normalizing the features, subtract mean divide by SD

pca_ = PCA(n_components=2)
#principalComponents_all_ = pca_.fit_transform(X_std)
#principalComponents_ = principalComponents_all_[:,:2]
principalComponents_ = pca_.fit_transform(X_std)

principalComponents_df = pd.DataFrame(principalComponents_, index=df.columns, columns = ['PC1', 'PC2'])

pc_df_reinxed = []
treatment_counts = []
for treatment in pt.treatments:
    pc_df_treatment = [x for x in principalComponents_df.index.to_list() if x[0]==treatment]
    pc_df_reinxed.extend(pc_df_treatment)

    treatment_counts.append(len(pc_df_treatment))

principalComponents_df = principalComponents_df.reindex(pc_df_reinxed)
principalComponents_np = principalComponents_df.values

# get dictionary of numpy arrays
taxon_pc_dict = {}
taxon_pc_labels = []
for taxon in taxa:
    taxon_index_names = [x for x in principalComponents_df.index.to_list() if x[1]==taxon]
    taxon_pc_df = principalComponents_df.loc[taxon_index_names]
    taxon_pc_dict[taxon] = taxon_pc_df.values
    taxon_pc_labels.extend(taxon_index_names)


treatment_counts = np.asarray(treatment_counts)

F = pt.get_F_2(principalComponents_np, treatment_counts)

# permute within each species
F_nested_list = []
for i in range(iter):
    #principalComponents_np_copy = np.copy(principalComponents_np)
    taxon_pc_dict_copy = taxon_pc_dict.copy()
    taxon_pc_permute_list = []
    for taxon in taxa:
        taxon_pc_dict_copy_i = taxon_pc_dict_copy[taxon]
        np.random.shuffle(taxon_pc_dict_copy_i)
        taxon_pc_permute_list.append(taxon_pc_dict_copy_i)

    principalComponents_np_permute = np.concatenate(taxon_pc_permute_list, axis=0)
    principalComponents_df_permute = pd.DataFrame(principalComponents_np_permute, index=taxon_pc_labels, columns = ['PC1', 'PC2'])
    # now resort by treatment
    pc_df_reinxed_permute = []
    for treatment in pt.treatments:
        pc_df_permute_treatment = [x for x in principalComponents_df_permute.index.to_list() if x[0]==treatment]
        pc_df_reinxed_permute.extend(pc_df_permute_treatment)

    principalComponents_df_permute = principalComponents_df_permute.reindex(pc_df_reinxed_permute)

    F_i = pt.get_F_2(principalComponents_df_permute.values, treatment_counts)
    F_nested_list.append(F_i)

p_value = (len([f for f in F_nested_list if f > F]) +1)/(iter+1)

sys.stderr.write("PERMANOVA F = %f, P = %f\n" % (round(F, 3), round(p_value, 6)))



# examine factor loadings
# variables were standardized, so we don't need to divive loadings by S.D. of variables
# if the variables weren't standardized, we'd end up w/ something like V * sqrt(E)/ std(X)
loadings = pca_.components_.T * np.sqrt(pca_.explained_variance_)

loading_squared_matrix = pd.DataFrame(loadings**2, columns=['PC1', 'PC2'], index=df.index)

loading_1_name=loading_squared_matrix.loc[loading_squared_matrix['PC1'].idxmax()].name
loading_1_rho_2=loading_squared_matrix.loc[loading_squared_matrix['PC1'].idxmax()][0]
loading_1_rho_2_idx = df.index.to_list().index(loading_1_name)

loading_2_name=loading_squared_matrix.loc[loading_squared_matrix['PC2'].idxmax()].name
loading_2_rho_2=loading_squared_matrix.loc[loading_squared_matrix['PC2'].idxmax()][1]
loading_2_rho_2_idx = df.index.to_list().index(loading_2_name)

transition_names = df.index.to_list()
# run permutation tests
permutation_dict = {}
for transition in transition_names:
    permutation_dict[transition] = {}
    permutation_dict[transition]['PC1'] = []
    permutation_dict[transition]['PC2'] = []

for i in range(iter_corr):
    X_copy = X.copy()
    #np.random.shuffle(pca_components_T_copy)
    for i in range(X_copy.shape[1]):
        np.random.shuffle(X_copy[:,i])

    X_copy_std=StandardScaler().fit_transform(X_copy) # normalizing the features, subtract mean divide by SD

    pca_copy_ = PCA(n_components=2)
    principalComponents_copy_ = pca_copy_.fit_transform(X_copy_std)
    loadings_copy_squared = (pca_copy_.components_.T * np.sqrt(pca_copy_.explained_variance_)) **2
    for loadings_row_idx, loadings_row in enumerate(loadings_copy_squared):

        permutation_dict[transition_names[loadings_row_idx]]['PC1'].append(loadings_row[0])
        permutation_dict[transition_names[loadings_row_idx]]['PC2'].append(loadings_row[1])


# get significant loadings
for index, row in loading_squared_matrix.iterrows():
    P_PC1 = (len([k for k in permutation_dict[transition]['PC1'] if k > float(row['PC1']) ])+1)/ (iter_corr+1)
    P_PC2 = (len([k for k in permutation_dict[transition]['PC2'] if k > float(row['PC2']) ])+1)/ (iter_corr+1)

    if P_PC1 < 0.05:
        sys.stderr.write("PC1: %s, rho^2 = %f, P=%f\n" % (index, round(row['PC1'], 3), round(P_PC1, 3)))

    if P_PC2 < 0.05:
        sys.stderr.write("PC2: %s, rho^2 = %f, P=%f\n" % (index, round(row['PC2'], 3), round(P_PC2, 3)))



# # now calculate pN/pS


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

for taxon in taxa:
    Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes(taxon)
    taxon_Lsyn_dict[taxon] = Lsyn
    taxon_Lnon_dict[taxon] = Lnon
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

        #print(synonymous_fmax_all)




fig = plt.figure(figsize = (9, 6))
gs = gridspec.GridSpec(nrows=2, ncols=3)
anova_pvalues = []
anova_F = []
mutation_spectra_list = [['AT_GC','AT_CG','AT_TA'],['GC_AT','GC_TA','GC_CG']]

axes = []

for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):

    for taxon_idx, taxon in enumerate(taxon_list):

        dnds_samples = []

        if taxon == 'J':
            treatments = ['0','2']
        else:
            treatments = pt.treatments

        #set_time = set_time_dict[taxon]

        #ax = fig.add_subplot(gs[taxon_idx*2:(taxon_idx*2)+2, taxon_list_idx])
        ax_pca = fig.add_subplot(gs[taxon_list_idx, taxon_idx ])
        axes.append(ax_pca)

        ax_pca.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )

        for treatment in treatments:

            PCs_ = principalComponents_df[principalComponents_df.index.str.contains(treatment+taxon)]

            ax_pca.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax_pca.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax_pca.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=2)

            ax_pca.scatter(PCs_.PC1.values, PCs_.PC2.values, \
                    c=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), s = 70, \
                    edgecolors=pt.get_colors(treatment), linewidth = 0.6, alpha = 0.8, zorder=4)#, edgecolors='none'

            pt.confidence_ellipse(PCs_.PC1.values, PCs_.PC2.values, ax_pca,
                n_std=2, edgecolor=pt.get_colors(treatment), linestyle='--', lw=4, zorder=3)

            # dn/ds
            populations_plot = [ treatment+taxon+replicate for replicate in replicates if treatment+taxon+replicate not in pt.populations_to_ignore ]
            taxon_treatment_dnds_appeared = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*taxon_Lsyn_dict[taxon]/taxon_Lnon_dict[taxon] for population in populations_plot]
            if len(taxon_treatment_dnds_appeared) < 2:
                continue
            dnds_samples.append(taxon_treatment_dnds_appeared)

        if taxon == 'J':
            fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1])

        else:
            fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1], dnds_samples[2])

        sys.stdout.write("%s: dN/dS one-way ANOVA: F = %g, P = %g\n" % (taxon, fvalue, pvalue))

        anova_pvalues.append(pvalue)
        anova_F.append(fvalue)

        ax_pca.text(-0.1, 1.07, sub_plot_labels[all_subplot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_pca.transAxes)

        all_subplot_counts += 1



#plt.plot([0.49, 0.49], [0, 1], color='k', ls=':', lw=5,transform=plt.gcf().transFigure, clip_on=False)


fig.text(0.5, -0.01, 'PC 1 (' + str(round(pca_.explained_variance_ratio_[0]*100,2)) + '%)', ha='center', va='center', fontsize=18)
fig.text(-0.01, 0.5, 'PC 2 (' + str(round(pca_.explained_variance_ratio_[1]*100,2)) + '%)', ha='center', va='center', rotation='vertical', fontsize=18)


custom_lines = [Line2D([0], [0], marker='o', markersize=10, color='w', markerfacecolor=pt.colors_dict['0']),
                Line2D([0], [0], marker='o', markersize=10, color='w', markerfacecolor=pt.colors_dict['1']),
                Line2D([0], [0], marker='o', markersize=10, color='w', markerfacecolor=pt.colors_dict['2'])]



axes[0].legend(custom_lines, ['1-day', '10-days', '100-day'], loc= 'upper right')


fig.subplots_adjust(hspace=0.4, wspace=0.6) #hspace=0.3, wspace=0.5
fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/mutation_spectra_pca.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()






reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(anova_pvalues, alpha=0.05, method='fdr_bh')

fig = plt.figure(figsize = (9, 6))
gs = gridspec.GridSpec(nrows=2, ncols=3)
all_subplot_counts = 0
dn_ds_count = 0
for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):
    for taxon_idx, taxon in enumerate(taxon_list):

        ax = fig.add_subplot(gs[taxon_list_idx, taxon_idx])
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold')
        dnds_samples = []
        for treatment in treatments:

            populations_plot = [ treatment+taxon+replicate for replicate in replicates if treatment+taxon+replicate not in pt.populations_to_ignore ]
            taxon_treatment_dnds_appeared = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*taxon_Lsyn_dict[taxon]/taxon_Lnon_dict[taxon] for population in populations_plot]
            if len(taxon_treatment_dnds_appeared) < 2:
                continue
            ax.scatter( [int(treatment)] * len(taxon_treatment_dnds_appeared), taxon_treatment_dnds_appeared,  marker=pt.plot_species_marker(taxon),  linewidth=2, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), s=100, zorder=2, alpha=0.8)
            if len(taxon_treatment_dnds_appeared) > 2:
                ax.errorbar(int(treatment),np.mean(taxon_treatment_dnds_appeared), yerr= 2*np.std(taxon_treatment_dnds_appeared) / np.sqrt(len(taxon_treatment_dnds_appeared)), linestyle='-', c = 'k', marker=pt.plot_species_marker(taxon), lw = 2.5,  zorder=3)
            #dnds_treatment.append(taxon_treatment_dnds_appeared)

            dnds_samples.append(taxon_treatment_dnds_appeared)

        ax.set_ylabel('pN/pS', fontsize = 12)

        ax.text(-0.1, 1.07, sub_plot_labels[all_subplot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)
        ax.text(0.7, 0.9, r'$F=$'+ str(round( anova_F[dn_ds_count],3) ), fontsize=10, ha='center', va='center', transform=ax.transAxes)
        ax.text(0.7, 0.8, r'$P_{BH}=$'+ str(round(pvals_corrected[dn_ds_count], 3)) , fontsize=10, ha='center', va='center', transform=ax.transAxes)

        all_subplot_counts+=1
        dn_ds_count += 1

        if taxon == 'J':
            ax.set_xticks([0,2])
            ax.set_xticklabels( ['1','100'] )
            ax.set_xlim([-0.3, 2.3])

            #fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1])
        else:
            ax.set_xticks([0,1,2])
            ax.set_xticklabels( ['1','10','100'] )
            ax.set_xlim([-0.3, 2.3])

            #fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1], dnds_samples[2])



fig.subplots_adjust(hspace=0.4, wspace=0.6) #hspace=0.3, wspace=0.5
fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/dn_ds.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()






#    #ma_mutation_spectrum = pt.get_ma_mutation_spectrum(taxon)
