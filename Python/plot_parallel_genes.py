from __future__ import division
import os, copy, sys, random
import matplotlib.pyplot as plt
#import matplotlib as mpl
from matplotlib import colors
from matplotlib.patches import Patch

import phylo_tools as pt
import parse_file

import scipy.stats as stats
import numpy as np


random.seed(123456789)

iter=1000

#taxa=['B', 'C', 'D', 'F', 'J', 'P']#pt.taxa
taxa=['C']
#taxa.remove('S')

treatments=pt.treatments
#treatments = ['0', '1', '2']

taxon_gene_dict = {}

for taxon in taxa:

    gene_dict = {}
    N_significant_genes = {}

    gene_data = parse_file.parse_gene_list(taxon)
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    locus_tag_to_gene_dict = {}
    for gene_name_idx, gene_name in enumerate(gene_names):
        gene = genes[gene_name_idx]
        if gene == '':
            continue
        locus_tag_to_gene_dict[gene_name] = genes[gene_name_idx]


    for treatment in treatments:

        genes_significant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_%ss_%s.txt" % ('gene', treatment+taxon))
        genes_nonsignificant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_not_significant_%ss_%s.txt" % ('gene', treatment+taxon))

        if os.path.exists(genes_significant_file_path) == False:
            continue

        genes_significant_file = open(genes_significant_file_path, 'r')
        first_line_significant = genes_significant_file.readline()

        N_significant_genes[treatment] = 0

        for line in genes_significant_file:
            line_split = line.strip().split(', ')
            gene_name = line_split[0]

            if gene_name not in gene_dict:
                gene_dict[gene_name] = {}

            gene_dict[gene_name][treatment] = 2

            N_significant_genes[treatment]+=1

        genes_significant_file.close()

        if os.path.exists(genes_nonsignificant_file_path) == True:
            genes_nonsignificant_file = open(genes_nonsignificant_file_path, 'r')
            first_line_nonsignificant = genes_nonsignificant_file.readline()

            for line in genes_nonsignificant_file:
                line_split = line.strip().split(', ')
                gene_name = line_split[0]

                if gene_name not in gene_dict:
                    gene_dict[gene_name] = {}

                gene_dict[gene_name][treatment] = 1

            genes_nonsignificant_file.close()

    # add zero for genes that you couldn't test
    for gene, gene_dict_i in gene_dict.items():
        for treatment in treatments:
            if treatment not in gene_dict_i:
                #if (treatment+taxon == '2J' ) or (treatment+taxon == '1F' ):
                #    gene_dict_i[treatment ] = -1
                #else:
                gene_dict_i[treatment ] = 0

    N_genes = len(gene_names)

    if len(N_significant_genes) == 3:

        intersect_1_10_100 = []
        intersect_1_10 = []
        intersect_1_100 = []
        intersect_10_100 = []

        for i in range(iter):

            # without replacement
            sample_1 = random.sample(set(range(N_genes)), N_significant_genes['0'])
            sample_10 = random.sample(set(range(N_genes)), N_significant_genes['1'])
            sample_100 = random.sample(set(range(N_genes)), N_significant_genes['2'])

            intersect_1_10_100.append(len(set(sample_1) & set(sample_10) & set(sample_100)))
            intersect_1_10.append(len(set(sample_1) & set(sample_10)))
            intersect_1_100.append(len(set(sample_1) & set(sample_100)))
            intersect_10_100.append(len(set(sample_10) & set(sample_100)))


        intersect_result_1_10 = len([key for key, value in gene_dict.items() if (gene_dict[key]['0']==2) and (gene_dict[key]['1']==2) ])
        intersect_result_1_100 = len([key for key, value in gene_dict.items() if (gene_dict[key]['0']==2) and (gene_dict[key]['2']==2) ])
        intersect_result_10_100 = len([key for key, value in gene_dict.items() if (gene_dict[key]['1']==2) and (gene_dict[key]['2']==2) ])
        intersect_result_1_10_100 = len([key for key, value in gene_dict.items() if (gene_dict[key]['0']==2) and (gene_dict[key]['1']==2) and (gene_dict[key]['2']==2)])

        # observed GREATER THAN expected

        intersect_result_1_10_p_value = (len([i for i in intersect_1_10 if i > intersect_result_1_10])+1) / (iter+1)
        intersect_result_1_100_p_value = (len([i for i in intersect_1_100 if i > intersect_result_1_100])+1) / (iter+1)
        intersect_result_10_100_p_value = (len([i for i in intersect_10_100 if i > intersect_result_10_100])+1) / (iter+1)
        intersect_result_1_10_100_p_value = (len([i for i in intersect_1_10_100 if i > intersect_result_1_10_100])+1) / (iter+1)


        sys.stderr.write("%s: 1-day & 10-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_10, intersect_result_1_10_p_value) )
        sys.stderr.write("%s: 1-day & 100-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_100, intersect_result_1_100_p_value) )
        sys.stderr.write("%s: 10-day & 100-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_10_100, intersect_result_10_100_p_value) )
        sys.stderr.write("%s: 1-day & 10-day & 100-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_10_100, intersect_result_1_10_100_p_value) )

    elif taxon == 'F':

        intersect_1_10 = []

        for i in range(iter):

            # without replacement
            sample_1 = random.sample(set(range(N_genes)), N_significant_genes['0'])
            sample_10 = random.sample(set(range(N_genes)), N_significant_genes['1'])
            intersect_1_10.append(len(set(sample_1) & set(sample_10)))

        intersect_result_1_10 = len([key for key, value in gene_dict.items() if (gene_dict[key]['0']==2) and (gene_dict[key]['1']==2) ])

        # observed GREATER THAN expected
        intersect_result_1_10_p_value = (len([i for i in intersect_1_10 if i > intersect_result_1_10])+1) / (iter+1)

        sys.stderr.write("%s: 1-day & 10-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_10, intersect_result_1_10_p_value) )


    elif taxon == 'J':

        intersect_1_100 = []

        for i in range(iter):

            # without replacement
            sample_1 = random.sample(set(range(N_genes)), N_significant_genes['0'])
            sample_100 = random.sample(set(range(N_genes)), N_significant_genes['2'])
            intersect_1_100.append(len(set(sample_1) & set(sample_100)))

        intersect_result_1_100 = len([key for key, value in gene_dict.items() if (gene_dict[key]['0']==2) and (gene_dict[key]['2']==2) ])

        # observed GREATER THAN expected
        intersect_result_1_100_p_value = (len([i for i in intersect_1_100 if i > intersect_result_1_100])+1) / (iter+1)

        sys.stderr.write("%s: 1-day & 10-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_100, intersect_result_1_100_p_value) )



    gene_dict_copy = copy.deepcopy(gene_dict)
    # remove genes if it's only significant in one treatment
    for gene in list(gene_dict_copy):
        if list(gene_dict_copy[gene].values()).count(2) < 2:
            del gene_dict_copy[gene]


    if taxon == 'J':
        for gene in list(gene_dict_copy):
            gene_dict_copy[gene]['1'] = -1

    #if taxon == 'F':
    #    for gene in list(gene_dict_copy):
    #        gene_dict_copy[gene]['2'] = -1

    gene_values = []
    gene_names = []
    for gene, gene_dict_i in sorted(gene_dict_copy.items()):

        #gene_dict_i_keys, gene_dict_i_values = sorted(gene_dict_i.items())
        gene_values.append(list( [gene_dict_i['0'], gene_dict_i['1'], gene_dict_i['2']]  ))
        if gene in locus_tag_to_gene_dict:
            gene_names.append(locus_tag_to_gene_dict[gene])
        else:
            gene_names.append(gene)

    taxon_gene_dict[taxon] = [gene_names,gene_values]



# 0 = no test
# 1 = tested, P > P*
# 2 = tested, P < P*

fig = plt.figure(figsize = (12, 24))
#fig = plt.figure(figsize=(6,3))
#fig, ax = plt.subplots(figsize=(2,7))

#cols = {0:'lightgrey',1:'orangered',2:'deepskyblue',-1:'black'}
cols = {0:'lightgrey',1:'orangered',2:'deepskyblue'}
cvr = colors.ColorConverter()
tmp = sorted(cols.keys())
cols_rgb = [cvr.to_rgb(cols[k]) for k in tmp]
intervals = np.asarray(tmp + [tmp[-1]+1]) - 0.5
cmap, norm = colors.from_levels_and_colors(intervals,cols_rgb)

legend_elements = [Patch(color='lightgrey', label=r'$n_{mut} < 3$'),
                    Patch(color='orangered', label=r'$P\nless P^{*}$'),
                    Patch(color='deepskyblue', label=r'$P < P^{*}$')]
#Patch(color='black', label="No data")

total_count = 0
for row_idx in range(2):

    for column_idx in range(3):

        ax = plt.subplot2grid((2, 3), (row_idx, column_idx), colspan=1)

        if total_count >= len(taxa):
            continue

        taxon = taxa[total_count]

        total_count+=1

        gene_names,gene_values = taxon_gene_dict[taxon]

        data = np.asarray(gene_values)

        ax.xaxis.tick_top()
        ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
        ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

        ax.set_yticklabels(gene_names, minor=False, fontsize=5)
        ax.set_xticklabels(['1-Day', '10-Days', '100-Days'], minor=False, fontsize=5.5)

        if total_count == 1:
            plt.legend(handles=legend_elements, bbox_to_anchor=(-0.65,1.12), loc="upper left", fontsize=8)

        plt.title(pt.latex_genus_bold_dict[taxon], fontsize=12)

        plt.pcolor(data,cmap = cmap, norm = norm, edgecolors='k', linewidths=1.5)

fig.subplots_adjust(wspace=.5, hspace=0.08)
fig.savefig(pt.get_path() + '/figs/gene_convergence_plot.pdf', format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
