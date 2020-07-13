from __future__ import division
import os, sys, copy, itertools, random
import numpy as np
from collections import Counter
import pandas as pd

import phylo_tools as pt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sklearn.decomposition import PCA

random.seed(123456789)

iter=1000

MCR = 1

treatments = ['0','1','2']
taxa = ['B', 'C', 'D', 'F', 'J', 'P']
maple_types = ['signature', 'complex', 'pathway', 'function']

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
                maple_annotation_dict[maple_name]['count'] = [taxon]
            else:
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

                    kegg_list.append(kegg_annotation)
                # get maple genes
                #    if kegg_to_maple_all[taxon][kegg_annotation] in maple_count_dict:
                #        maple_count_dict[kegg_to_maple_all[taxon][kegg_annotation]] += 1
                    maple_list.append(kegg_to_maple_all[taxon][kegg_annotation])
                #    else:
                #        maple_count_dict[kegg_to_maple_all[taxon][kegg_annotation]] = 1

        significant_genes.close()

        kegg_taxon_treatment[treatment+taxon] = kegg_list
        maple_taxon_treatment[treatment+taxon] = list(set(maple_list))



# randomly sample from the set of genes that CAN be mapped
null_intersection_size_dict = {}

for treatment in treatments:

    null_intersection_size_dict[treatment] = []

    if (treatment == '1'):
        taxa_ =  ['B','C','D','F','P']
    elif (treatment == '2'):
        taxa_ =  ['B','C','D','J','P']
    else:
        taxa_ =  pt.taxa

    #for i in range(1, len(taxa_)+1):
    #    null_intersection_size_dict[treatment][i] = []

    for i in range(iter):
        null_maple_dict = {}

        for taxon_ in taxa_:

            N = len(kegg_taxon_treatment[treatment+taxon_])
            sample_kegg = random.sample(list(kegg_to_maple_all[taxon_].keys()), N)
            sample_maple = list(set([kegg_to_maple_all[taxon_][kegg_] for kegg_ in sample_kegg]))
            null_maple_dict[taxon_] = sample_maple

        null_decay_curve = []
        for j in range(1, len(taxa_)+1):

            comb_taxa = list(itertools.combinations(taxa_, j))
            sum_j = 0

            for comb_taxa_j in comb_taxa:

                all_maple_modules = []
                for comb_taxa_j_k in comb_taxa_j:
                    all_maple_modules.extend(null_maple_dict[comb_taxa_j_k])
                maple_count = Counter(all_maple_modules)

                intersection_size = list(dict(maple_count).values()).count(j)

                sum_j += intersection_size

            #null_intersection_size_dict[treatment][j].append(sum_j)

            null_decay_curve.append(sum_j)

        null_intersection_size_dict[treatment].append(null_decay_curve)



observed_intersection_dict = {}
# get the observed decay curve
# plot each of them for each treatment w/ 10^4 null simulations
for treatment in treatments:

    observed_intersection_dict[treatment] = []

    if (treatment == '1'):
        taxa_ =  ['B','C','D','F','P']
    elif (treatment == '2'):
        taxa_ =  ['B','C','D','J','P']
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


maple_matrix = pd.DataFrame.from_dict(maple_taxon_treatment_dict)
maple_matrix = maple_matrix.fillna(0)
maple_matrix = maple_matrix.T

maple_matrix_centered = maple_matrix - np.mean(maple_matrix, axis = 0)

pca = PCA(n_components = 2)
X_pca = pca.fit_transform(maple_matrix_centered)

#pca = PCA(n_components = 2).fit(maple_matrix_centered)
#df_out = pca.fit_transform(X)
df_pc = pd.DataFrame(data=X_pca, index=maple_matrix.index.values)


fig = plt.figure(figsize = (10, 6))
gs = gridspec.GridSpec(nrows=3, ncols=4)
sub_plot_labels = ['a', 'b', 'c']
for treatment_idx, treatment in enumerate(treatments):

    #ax = plt.subplot2grid((3, 5), (treatment_idx, 0), colspan=2, rowspan=1)
    ax = fig.add_subplot(gs[treatment_idx, 0])

    for null_i in null_intersection_size_dict[treatment]:
        ax.plot(range(1, len(null_i)+1), null_i, '-', alpha=0.1, zorder=1, c = pt.get_colors(treatment))

    ax.plot(range(1, len(observed_intersection_dict[treatment])+1), observed_intersection_dict[treatment], 'o-', alpha=1, c='k',lw=2,zorder=2)

    #ax.set_yscale('log', basey=10)
    ax.set_xlim(0.3, 6.4)
    #ax.set_title( str(10**int(treatment)) + '-day', fontsize=20 )
    if treatment == '0':
        title = '1-Day'
    elif treatment == '1':
        title = '10-Days'
    elif treatment == '2':
        title = '100-Days'
        ax.set_xlabel('Number of taxa', fontsize = 16)
    else:
        continue
    ax.text(0.4, 0.8, title, fontsize=12, transform=ax.transAxes)

    ax.text(0, 1.1, sub_plot_labels[treatment_idx], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)


# now plot PCA

#ax_pca = plt.subplot2grid((3, 5), (0, 1), colspan=3, rowspan=3)
ax_pca = fig.add_subplot(gs[0:3, 1:4])

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



ax_pca.set_xlabel('PC 1 (' + str(round(pca.explained_variance_ratio_[0]*100,2)) + '%)' , fontsize = 16)
ax_pca.set_ylabel('PC 2 (' + str(round(pca.explained_variance_ratio_[1]*100,2)) + '%)' , fontsize = 16)

#### perform PERMANOVA
F, p = pt.run_permanova(df_pc.values, number_treatment_reps)
ax_pca.text(0.8, 0.9, r'$F$=' + str(round(F,2)), fontsize = 14, transform=ax_pca.transAxes)
ax_pca.text(0.8, 0.83, r'$P\nless 0.05$', fontsize = 14, transform=ax_pca.transAxes)

ax_pca.text(0, 1.1, 'd', fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_pca.transAxes)

fig.text(0.05, 0.5, 'Number of MAPLE modules with\nsignificant multiplicity', ha='center', va='center', rotation='vertical', fontsize=15)

fig_name = pt.get_path() + '/figs/convergence_decay.pdf'
fig.subplots_adjust(wspace=0.55)
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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
        out_2 = str(len(maple_i_dict['2'])) + '/5'
    else:
        out_2 = '0/5'

    convergence_table.write(", ".join([maple_i, type, description, out_0, out_1, out_2 ]) + '\n' )

convergence_table.close()
