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

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

np.random.seed(123456789)

iter=1

spectra_latex_dict = {'AT_GC':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{G}:\mathrm{C}$',
                    'AT_CG':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{C}:\mathrm{G}$',
                    'AT_TA':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_AT':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{A}:\mathrm{T}$',
                    'GC_TA':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_CG':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{C}:\mathrm{G}$'}


treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']
iter
spectra_dict = {}

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

X=StandardScaler().fit_transform(X) # normalizing the features, subtract mean divide by SD

normalised_df = pd.DataFrame(X,columns=df.index, index=df.columns)
pca_ = PCA(n_components=2)
principalComponents_ = pca_.fit_transform(X)

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


sys.stderr.write("PERMANOVA F = %f, P = %f\n" % (round(F, 3), round(p_value, 3)))




#fig, ax = plt.subplots(figsize=(8, 6))

fig = plt.figure(figsize = (12, 9))
gs = gridspec.GridSpec(nrows=6, ncols=4)


sub_plot_counts = 0

for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):

    for taxon_idx, taxon in enumerate(taxon_list):

        if taxon == 'J':
            treatments = ['0','2']
        else:
            treatments = pt.treatments

        #set_time = set_time_dict[taxon]

        ax = fig.add_subplot(gs[taxon_idx*2:(taxon_idx*2)+2, 2+taxon_list_idx])
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )

        for treatment in treatments:

            PCs_ = principalComponents_df[principalComponents_df.index.str.contains(treatment+taxon)]

            ax.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=2)

            ax.scatter(PCs_.PC1.values, PCs_.PC2.values, \
                    c=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), s = 70, \
                    edgecolors=pt.get_colors(treatment), linewidth = 0.6, alpha = 0.8, zorder=3)#, edgecolors='none'

            pt.confidence_ellipse(PCs_.PC1.values, PCs_.PC2.values, ax,
                n_std=2, edgecolor=pt.get_colors(treatment), linestyle='--', lw=4)

# mutation spectrum
#ax_snp_1 = fig.add_subplot(gs[0, 0:2])
#ax_snp_2 = fig.add_subplot(gs[1, 0:2])

mutation_spectra_list = [['AT_GC','AT_CG','AT_TA'],['GC_AT','GC_TA','GC_CG']]
#ax_snp_1.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
#ax_snp_2.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)

count = 0
for snp_i in range(2):
    ax_snp = fig.add_subplot(gs[snp_i*2:(snp_i+1)*2, 0:2])
    ax_snp.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
    ax_snp.set_xticklabels([])
    ax_snp.set_xticks([])
    for spectrum_idx, spectrum in enumerate(mutation_spectra_list[snp_i]):
        for taxon in ['B','C','D','J']:
            for treatment in pt.treatments:
                spectrum_reps = []
                for replicate in replicates:
                    population = treatment+taxon+replicate
                    if population in spectra_dict:
                        spectrum_reps.append(spectra_dict[population][spectrum])

                if len(spectrum_reps) < 2:
                    continue

                spectrum_reps = np.asarray(spectrum_reps) - pt.get_ma_mutation_spectrum(taxon)[spectrum]

                ax_snp.scatter(list(itertools.repeat(count, len(spectrum_reps))), spectrum_reps,\
                        c=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), s = 30, \
                        edgecolors=pt.get_colors(treatment), linewidth = 0.6, alpha = 0.8, zorder=2)#, edgecolors='none'


                count+=1


        ax_snp.text((spectrum_idx+0.5)*(1/3), -0.05, spectra_latex_dict[spectrum], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_snp.transAxes)


        if spectrum_idx < 2:

            ax_snp.axvline(x=count-0.5, color='k', linestyle='--', alpha = 0.8, zorder=1)




fig.text(0.75, -0.01, 'PC 1 (' + str(round(pca_.explained_variance_ratio_[0]*100,2)) + '%)', ha='center', va='center', fontsize=18)
fig.text(0.5, 0.5, 'PC 2 (' + str(round(pca_.explained_variance_ratio_[1]*100,2)) + '%)', ha='center', va='center', rotation='vertical', fontsize=18)

fig.text(-0.01, 0.5, 'Observed spectra - MA spectra', ha='center', va='center', rotation='vertical', fontsize=18)


fig.subplots_adjust(hspace=0.4, wspace=0.38) #hspace=0.3, wspace=0.5
fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/mutation_spectra_pca.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



#    #ma_mutation_spectrum = pt.get_ma_mutation_spectrum(taxon)
