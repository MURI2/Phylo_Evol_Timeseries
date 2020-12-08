from __future__ import division
import os, sys, copy, itertools, random, math
import numpy as np
from collections import Counter
from itertools import combinations

import scipy.stats as stats
import pandas as pd

import phylo_tools as pt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from sklearn.decomposition import PCA


import ete3

random.seed(123456789)
np.random.seed(123456789)

significant_digits=3

iter=10000
#iter=100

regression_permutations = 100000
#regression_permutations = 100

MCR = 1


treatments = ['0','1','2']
taxa = ['B', 'C', 'D', 'F', 'J', 'P']
maple_types = ['signature', 'complex', 'pathway', 'function']


# Get phylogenetic distance

tree = ete3.Tree(pt.get_path() + '/data/tree/RAxML_bipartitions.Task2', quoted_node_names=True, format=1)

phylogenetic_distance_dict = {}



for taxa_pair in combinations(taxa, 2):

    phylogenetic_distance =  tree.get_distance(pt.tree_name_dict[taxa_pair[0]] , pt.tree_name_dict[taxa_pair[1]] )
    phylogenetic_distance_dict[taxa_pair] = phylogenetic_distance


# Now work on kegg genes

maple_annotation_dict = {}

kegg_to_maple_all = {}

maple_to_kegg_all = {}
protein_to_kegg_all = {}
# get list of maple modules in all taxa
for taxon in taxa:
    # map KEGG genes onto pathays
    # get pathways
    # go through signatures, complexes, pathways, and functions
    kegg_to_maple_all[taxon] = {}
    maple_to_kegg_all[taxon] = {}
    protein_to_kegg_all[taxon] = {}
    protein_id_kegg = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/query.fst.ko' % taxon, 'r')
    # make protein ID => KEGG map
    for i, line in enumerate(protein_id_kegg):
        line = line.strip()
        items = line.split("\t")
        protein_id = items[0]
        if items[1] != 'K_NA':
            protein_to_kegg_all[taxon][items[0]] = items[1]

    for maple_type in maple_types:
        maple_file = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/module_%s.tsv' % (taxon, maple_type) , 'r')
        first_line_maple = maple_file.readline()
        first_line_maple = first_line_maple.strip()
        first_line_items_maple = first_line_maple.split("\t")
        for i, line in enumerate(maple_file):
            line = line.strip()
            items = line.split("\t")

            if float(items[10]) < MCR:
                continue
            # remove rows that are less than 80% complete
            # query(coverage) = MCR % (ITR)
            # query(coverage/max) = MCR % (WC)
            # query(coverage/mode) = Q-value
            maple_name = items[0].split('_')[0]
            #maple_to_keep[maple_name] = {}#[items[1], items[2]]
            #maple_to_keep[maple_name]['type'] = items[1]
            #maple_to_keep[maple_name]['description'] = items[2]
            #      #kegg_to_maple_all[taxon][matrix_items[0]] =

            if maple_name not in maple_annotation_dict:
                maple_annotation_dict[maple_name] = {}
                maple_annotation_dict[maple_name]['type'] = items[1]
                maple_annotation_dict[maple_name]['description'] = items[2]
                maple_annotation_dict[maple_name]['count'] = []

            if taxon not in maple_annotation_dict[maple_name]['count']:

                maple_annotation_dict[maple_name]['count'].append(taxon)

            maple_matrix = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/KAAS/%s_matrix.txt' % (taxon, maple_name) , 'r')
            maple_kegg = []
            # only get kegg genes that exist in the taxon's genome
            for matrix_line_i, matrix_line in enumerate(maple_matrix):
                matrix_line = matrix_line.strip()
                matrix_items = matrix_line.split('\t')
                if 'K' not in matrix_items[0]:
                    continue
                if len(matrix_items) < 3:
                    continue
                if int(matrix_items[1]) == 0:
                    continue
                maple_kegg.append(matrix_items[0])

            maple_to_kegg_all[taxon][maple_name] = maple_kegg

    # now map kegg to pathways
    #kegg_to_maple_for_simulation[treatment+taxon] = []
    #for kegg_i in kegg_list:

    for maple, kegg_genes in maple_to_kegg_all[taxon].items():
        for kegg_gene in kegg_genes:
            kegg_to_maple_all[taxon][kegg_gene] = maple


# filter kegg_to_maple_all to only include MAPLE modules
# with maple_annotation_dict[maple_name]['count'] = 6




kegg_taxon_treatment = {}
maple_taxon_treatment = {}
for treatment in treatments:
    # make protein ID => KEGG map
    for taxon in taxa:
        significant_genes_path=pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significant_genes_path) == False:
            continue
        significant_genes = open(significant_genes_path, 'r')

        first_line = significant_genes.readline()
        first_line = first_line.strip()
        first_line_items = first_line.split(",")

        kegg_list = []
        maple_list = []
        # get KEGG annotation of genes with significant multiplicity
        for i, line in enumerate(significant_genes):
            line = line.strip()
            items = line.split(",")
            protein_id = items[1].strip()
            if protein_id in protein_to_kegg_all[taxon]:
                kegg_annotation = protein_to_kegg_all[taxon][protein_id]
                if kegg_annotation in kegg_to_maple_all[taxon]:

                    maple_annotation = kegg_to_maple_all[taxon][kegg_annotation]
                    if len(maple_annotation_dict[maple_annotation]['count']) == 6:

                        kegg_list.append(kegg_annotation)

                        maple_list.append(maple_annotation)

        significant_genes.close()

        kegg_taxon_treatment[treatment+taxon] = list(set(kegg_list))
        maple_taxon_treatment[treatment+taxon] = list(set(maple_list))



# randomly sample from the set of genes that CAN be mapped
null_intersection_size_maple_dict = {}
null_intersection_size_kegg_dict = {}

for treatment in treatments:

    null_intersection_size_maple_dict[treatment] = []
    null_intersection_size_kegg_dict[treatment] = []

    if (treatment == '1'):
        taxa_ =  ['B','C','D','F','P']
    #elif (treatment == '2'):
    #    taxa_ =  ['B','C','D','J','P']
    else:
        taxa_ =  pt.taxa

    for i in range(iter):
        null_kegg_dict = {}
        null_maple_dict = {}

        for taxon_ in taxa_:

            N = len(kegg_taxon_treatment[treatment+taxon_])
            sample_kegg = random.sample(list(kegg_to_maple_all[taxon_].keys()), N)
            sample_maple = list(set([kegg_to_maple_all[taxon_][kegg_] for kegg_ in sample_kegg]))
            null_kegg_dict[taxon_] = sample_kegg
            null_maple_dict[taxon_] = sample_maple

        null_decay_curve_kegg = []
        null_decay_curve_maple = []
        for j in range(1, len(taxa_)+1):

            comb_taxa = list(itertools.combinations(taxa_, j))
            sum_j_kegg = 0
            sum_j_maple = 0

            for comb_taxa_j in comb_taxa:

                all_kegg_modules = []
                all_maple_modules = []
                for comb_taxa_j_k in comb_taxa_j:
                    all_kegg_modules.extend(null_kegg_dict[comb_taxa_j_k])
                    all_maple_modules.extend(null_maple_dict[comb_taxa_j_k])

                kegg_count = Counter(all_kegg_modules)
                maple_count = Counter(all_maple_modules)

                intersection_size_kegg = list(dict(kegg_count).values()).count(j)
                intersection_size_maple = list(dict(maple_count).values()).count(j)

                sum_j_kegg += intersection_size_kegg
                sum_j_maple += intersection_size_maple

            #null_intersection_size_dict[treatment][j].append(sum_j)
            null_decay_curve_kegg.append(sum_j_kegg)
            null_decay_curve_maple.append(sum_j_maple)

        null_intersection_size_kegg_dict[treatment].append(null_decay_curve_kegg)
        null_intersection_size_maple_dict[treatment].append(null_decay_curve_maple)


observed_intersection_dict = {}
# get the observed decay curve
# plot each of them for each treatment w/ 10^4 null simulations
for treatment in treatments:

    observed_intersection_dict[treatment] = []

    if (treatment == '1'):
        taxa_ =  ['B','C','D','F','P']
    else:
        taxa_ =  pt.taxa

    for j in range(1, len(taxa_)+1):
        #for taxon_ in taxa_:

        comb_taxa = list(itertools.combinations(taxa_, j))
        sum_j = 0

        for comb_taxa_j in comb_taxa:

            all_maple_modules = []
            for comb_taxa_j_k in comb_taxa_j:
                all_maple_modules.extend(maple_taxon_treatment[treatment+comb_taxa_j_k])
            maple_count = Counter(all_maple_modules)

            intersection_size = list(dict(maple_count).values()).count(j)

            sum_j += intersection_size

        observed_intersection_dict[treatment].append(sum_j)


maple_taxon_treatment_dict = {}
for key, value in maple_taxon_treatment.items():
    maple_taxon_treatment_dict[key] = {}
    for value_i in value:
        maple_taxon_treatment_dict[key][value_i] = 1


kegg_taxon_treatment_dict = {}
for key, value in kegg_taxon_treatment.items():
    kegg_taxon_treatment_dict[key] = {}
    for value_i in value:
        kegg_taxon_treatment_dict[key][value_i] = 1



maple_matrix = pd.DataFrame.from_dict(maple_taxon_treatment_dict)
maple_matrix = maple_matrix.fillna(0)
maple_matrix = maple_matrix.T

maple_matrix_centered = maple_matrix - np.mean(maple_matrix, axis = 0)

pca = PCA(n_components = 2)
X_pca = pca.fit_transform(maple_matrix_centered)

#pca = PCA(n_components = 2).fit(maple_matrix_centered)
#df_out = pca.fit_transform(X)
df_pc = pd.DataFrame(data=X_pca, index=maple_matrix.index.values)


fig = plt.figure(figsize = (7, 6))
#gs = gridspec.GridSpec(nrows=3, ncols=4)
gs = gridspec.GridSpec(nrows=6, ncols=5)

sub_plot_labels = ['a', 'b', 'c']
for treatment_idx, treatment in enumerate(treatments):

    #ax = plt.subplot2grid((3, 5), (treatment_idx, 0), colspan=2, rowspan=1)

    ax = fig.add_subplot(gs[treatment_idx*2:(treatment_idx*2)+2  , 0:2])

    #ax = fig.add_subplot(gs[treatment_idx, 0])

    for null_i in null_intersection_size_maple_dict[treatment]:
        ax.plot(range(1, len(null_i)+1), null_i, '-', alpha=0.1, zorder=1, c = pt.get_colors(treatment))

    ax.plot(range(1, len(observed_intersection_dict[treatment])+1), observed_intersection_dict[treatment], 'o-', alpha=1, c='k',lw=2,zorder=2)

    #ax.set_yscale('log', basey=10)
    ax.set_xlim(0.8, 6.4)
    #ax.set_title( str(10**int(treatment)) + '-day', fontsize=20 )
    if treatment == '0':
        title = '1-Day'
    elif treatment == '1':
        title = '10-Days'
    elif treatment == '2':
        title = '100-Days'
        ax.set_xlabel('Number of taxa', fontsize = 14)
    else:
        continue
    #ax.text(0.62, 0.85, title, fontsize=10, transform=ax.transAxes)

    ax.text(0, 1.1, sub_plot_labels[treatment_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax.transAxes)



    #legend_elements = [Line2D([0], [0], color='k', ls='-', lw=1.5, label='Observed'),
    #                Line2D([0], [0], color=pt.get_colors(treatment), ls='-', lw=1.5, label='Null')]
    legend_elements = [Line2D([0], [0], marker='o', color='k', markerfacecolor='k', markersize=6, label='Observed'),
                    Line2D([0], [0], color=pt.get_colors(treatment), ls='-', lw=1.5, label='Null')]

    ax.legend(handles=legend_elements, loc='upper right',  prop={'size': 8})

    ax.set_title(title, fontdict={'fontsize': 10 })
    # 'fontweight': 'medium'


# now plot PCA


ax_phylo = fig.add_subplot(gs[0:3, 2:])

ax_phylo.text(0, 1.05, 'd', fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)

treatment_distance_dict = {}


treatments_all = []
relative_intersection_size_all = []

for treatment_idx, treatment in enumerate(treatments):

    if treatment == '1':
        taxa_iter = taxa
        taxa_iter.remove('J')

    else:
        taxa_iter = taxa

    distances = []
    relative_intersection_size = []

    for taxa_pair in combinations(taxa_iter, 2):

        #maple_1 = set(maple_taxon_treatment[treatment+taxa_pair[0]])
        #maple_2 = set(maple_taxon_treatment[treatment+taxa_pair[1]])
        maple_1 = set(maple_taxon_treatment[treatment+taxa_pair[0]])
        maple_2 = set(maple_taxon_treatment[treatment+taxa_pair[1]])

        distances.append(tree.get_distance(pt.tree_name_dict[taxa_pair[0]] , pt.tree_name_dict[taxa_pair[1]]) )

        relative_intersection_size.append(len(maple_1.intersection(maple_2)  ) / len(maple_1.union(maple_2)))


    treatments_all += distances
    relative_intersection_size_all += relative_intersection_size


    distances = np.asarray(distances)
    relative_intersection_size = np.asarray(relative_intersection_size)

    ax_phylo.scatter(distances, relative_intersection_size, c = pt.get_colors(treatment), s=50, alpha=0.8)

    #relative_intersection_size_nonzero_idx = relative_intersection_size>0
    #distances_nonzero = distances[relative_intersection_size_nonzero_idx]
    #relative_intersection_size_nonzero = relative_intersection_size[relative_intersection_size_nonzero_idx]

    slope, intercept, r_value, p_value, std_err = stats.linregress(distances, relative_intersection_size)
    slope_permuted_list = []
    for permute in range(regression_permutations):
        distances_permuted = np.random.permutation(distances)
        relative_intersection_size_permuted = np.random.permutation(relative_intersection_size)
        slope_permuted, intercept_permuted, r_value_permuted, p_valu_permutede, std_err_permuted = stats.linregress(distances_permuted, relative_intersection_size_permuted)
        slope_permuted_list.append(slope_permuted)

    slope_permuted_list = np.asarray(slope_permuted_list)

    x_range =  np.linspace(0.3, 0.7, 10000)
    y_fit_range = (slope*x_range + intercept)

    p_value_permuted = len(slope_permuted_list[slope_permuted_list<slope]) / regression_permutations

    sys.stdout.write("%d-day distance-decay: slope=%f, P=%f\n" % (10**int(treatment), slope, p_value_permuted))

    y_position_labels = [0.9,0.73,0.56]

    position_dict_x = {'0':0.71, '1':0.73, '2':0.74}
    position_dict_y = {'0':0.92, '1':0.80, '2':0.68}

    if p_value_permuted < 0.05:
        ax_phylo.plot(x_range, y_fit_range, c=pt.get_colors(treatment), lw=2.5, linestyle='--', zorder=2)

    #    #ax_phylo.text(0.8,0.9, r'$y \sim x^{{{}}}$'.format(str( round(slope, 3) )), fontsize=12, color=pt.get_colors(treatment), ha='center', va='center', transform=ax_phylo.transAxes  )
    #    ax_phylo.text(position_dict_x[treatment], position_dict_y[treatment], r'$y \sim x^{{{}}} \; P <0.05$'.format(str( round(slope, 3) )), fontsize=11, color=pt.get_colors(treatment), ha='center', va='center', transform=ax_phylo.transAxes  )
    #else:
    #    ax_phylo.text(position_dict_x[treatment], position_dict_y[treatment], r'$y \sim x^{{{}}} \; P \nless 0.05$'.format(str( round(slope, 3) )), fontsize=11, color=pt.get_colors(treatment), ha='center', va='center', transform=ax_phylo.transAxes  )
    rounded_p_value_permuted =  round(p_value_permuted, significant_digits - int(math.floor(math.log10(abs(p_value_permuted)))) - 1)

    ax_phylo.text(position_dict_x[treatment], position_dict_y[treatment], r'$y \sim x^{{{}}} \; P = {{{}}}$'.format(str( round(slope, 3) ),  str(rounded_p_value_permuted) ) , fontsize=11, color=pt.get_colors(treatment), ha='center', va='center', transform=ax_phylo.transAxes  )



    treatment_distance_dict[treatment] = {}
    treatment_distance_dict[treatment]['distances'] = distances
    treatment_distance_dict[treatment]['relative_intersection_size'] = relative_intersection_size




ax_phylo.set_xlabel('Phylogenetic distance' , fontsize = 10)
ax_phylo.set_ylabel('Jaccard index of MAPLE modules' , fontsize = 10)

treatments_all = np.asarray(treatments_all)
relative_intersection_size_all = np.asarray(relative_intersection_size_all)


slope_x, intercept_x, r_value_x, p_value_x, std_err_x = stats.linregress(treatments_all, relative_intersection_size_all)
slopes_nullll = []
for permute_x in range(regression_permutations):
    treatments_all_permute = np.random.permutation(treatments_all)
    relative_intersection_size_all_permute = np.random.permutation(relative_intersection_size_all)

    slope_xx, intercept_xx, r_value_xx, p_value_xx, std_err_xx = stats.linregress(treatments_all_permute, relative_intersection_size_all_permute)

    slopes_nullll.append(slope_xx)

slopes_nullll = np.asarray(slopes_nullll)



# slope difference test
for treatment_pair in combinations(treatments, 2):

    distances_1 = treatment_distance_dict[treatment_pair[0]]['distances']
    distances_2 = treatment_distance_dict[treatment_pair[1]]['distances']

    relative_intersection_size_1 = treatment_distance_dict[treatment_pair[0]]['relative_intersection_size']
    relative_intersection_size_2 = treatment_distance_dict[treatment_pair[1]]['relative_intersection_size']

    merged_variables = list(zip(distances_1,relative_intersection_size_1)) + list(zip(distances_2,relative_intersection_size_2))
    merged_variables = np.asarray(merged_variables)

    sample_size_1 = len(distances_1)
    sample_size_2 = len(distances_2)

    slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(distances_1, relative_intersection_size_1)
    slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(distances_2, relative_intersection_size_2)


    slope_difference = abs(slope_1 - slope_2)

    slope_differences_null = []
    for permute in range(regression_permutations):

        merged_variables_permuted = np.random.permutation(merged_variables)

        merged_variables_permuted_1 = merged_variables_permuted[:sample_size_1]
        merged_variables_permuted_2 = merged_variables_permuted[sample_size_1+1:]

        distances_permuted_1 = [l[0] for l in merged_variables_permuted_1]
        relative_intersection_size_permuted_1 = [l[1] for l in merged_variables_permuted_1]

        distances_permuted_2 = [l[0] for l in merged_variables_permuted_2]
        relative_intersection_size_permuted_2 = [l[1] for l in merged_variables_permuted_2]

        slope_permute_1, intercept_permute_1, r_value_permute_1, p_value_permute_1, std_err_permute_1 = stats.linregress(distances_permuted_1, relative_intersection_size_permuted_1)

        slope_permute_2, intercept_permute_2, r_value_permute_2, p_value_permute_2, std_err_permute_2 = stats.linregress(distances_permuted_2, relative_intersection_size_permuted_2)

        slope_differences_null.append(abs(slope_permute_1-slope_permute_2))

    slope_differences_null = np.asarray(slope_differences_null)

    P =  len(slope_differences_null[slope_differences_null>slope_difference]) / regression_permutations

    sys.stdout.write("%d vs %d-day slope difference test: |slope differnce|=%f, P=%f\n" % (10**int(treatment_pair[0]), 10**int(treatment_pair[1]), slope_difference, P))









#ax_pca = plt.subplot2grid((3, 5), (0, 1), colspan=3, rowspan=3)
ax_pca = fig.add_subplot(gs[3:6, 2:])

ax_pca.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
ax_pca.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
ax_pca.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=2)

for index, row in df_pc.iterrows():

    ax_pca.scatter(row[0], row[1], marker=pt.plot_species_marker(index[1]), s = 220, \
        linewidth=3, facecolors=pt.get_scatter_facecolor(index[1], index[0]), \
        edgecolors='none', alpha=0.8, zorder=3)

# plot ellipses
number_treatment_reps = []
for treatment in treatments:

    treatment_df = df_pc[df_pc.index.str.startswith(treatment) ]

    number_treatment_reps.append(treatment_df.shape[0])

    pt.confidence_ellipse(treatment_df.values[:,0],treatment_df.values[:,1], ax_pca,
        n_std=2, edgecolor=pt.get_colors(treatment), linestyle='--', lw=3)



ax_pca.set_xlabel('PC 1 (' + str(round(pca.explained_variance_ratio_[0]*100,2)) + '%)' , fontsize = 12)
ax_pca.set_ylabel('PC 2 (' + str(round(pca.explained_variance_ratio_[1]*100,2)) + '%)' , fontsize = 12)

#### perform PERMANOVA
F, P = pt.run_permanova(df_pc.values, number_treatment_reps)
rounded_F =  round(F, significant_digits - int(math.floor(math.log10(abs(F)))) - 1)
rounded_P =  round(P, significant_digits - int(math.floor(math.log10(abs(P)))) - 1)

ax_pca.text(0.77, 0.9, r'$F$=' + str(round(rounded_F,3)), fontsize = 11, transform=ax_pca.transAxes)
#ax_pca.text(0.8, 0.8, r'$P\nless 0.05$', fontsize = 11, transform=ax_pca.transAxes)
ax_pca.text(0.765, 0.8, r'$P=$' +str(round(rounded_P,3)), fontsize = 11, transform=ax_pca.transAxes)

ax_pca.text(0, 1.05, 'e', fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_pca.transAxes)

fig.text(0.05, 0.5, 'Number of intersecting MAPLE modules\nwith significant multiplicity', ha='center', va='center', rotation='vertical', fontsize=14)

fig_name = pt.get_path() + '/figs/convergence_decay.pdf'
fig.subplots_adjust(wspace=1, hspace=1.2)
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()




# reformat data for output
maple_dict_count = {}
for key, maple_list in maple_taxon_treatment.items():

    for maple_i in maple_list:

        if maple_i not in maple_dict_count:
            maple_dict_count[maple_i] = {}

        if key[0] not in maple_dict_count[maple_i]:
            maple_dict_count[maple_i][key[0]] = []

        if key[1] not in maple_dict_count[maple_i][key[0]]:
            maple_dict_count[maple_i][key[0]].append(key[1])


# print it out to a file
# idk reformat and do it later

convergence_table = open(pt.get_path() + '/data/parallel_pathways.txt', 'w')
# alot of ribosome complex, merge those first
convergence_table.write(", ".join(["Module ID", "Module type", "Description", "1-day", "10-day", "100-day" ]) + '\n' )

for maple_i, maple_i_dict in maple_dict_count.items():
    keep_maple=False
    for treatment in treatments:
        if treatment in  maple_i_dict:
            if len(maple_i_dict[treatment]) >= 2:
                keep_maple=True

    if keep_maple == False:
        continue

    type = maple_annotation_dict[maple_i]['type']
    description = maple_annotation_dict[maple_i]['description']
    description = description.replace(',', ';')

    description = description.split(';')[0]
    description = description.replace(' (Pentose phosphate cycle)', '')
    description = description.replace(' (PE)', '')

    #if 'Ribosome' in description:
    #    continue

    if '0' in maple_i_dict:
        out_0 = str(len(maple_i_dict['0'])) + '/6'
    else:
        out_0 = '0/6'

    if '1' in maple_i_dict:
        out_1 = str(len(maple_i_dict['1'])) + '/5'
    else:
        out_1 = '0/5'

    if '2' in maple_i_dict:
        out_2 = str(len(maple_i_dict['2'])) + '/6'
    else:
        out_2 = '0/6'

    convergence_table.write(", ".join([maple_i, type, description, out_0, out_1, out_2 ]) + '\n' )

convergence_table.close()
