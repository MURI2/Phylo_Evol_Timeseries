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

taxa=pt.taxa

treatments=pt.treatments
#treatments = ['0', '1', '2']

taxon_gene_dict = {}

for taxon in taxa:

    if taxon == 'J':
        continue

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

        sys.stderr.write("%s: 1-day & 100-day %d genes, p = %f\n" % (pt.genus_dict[taxon] ,intersect_result_1_100, intersect_result_1_100_p_value) )



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

cols = {0:'lightgrey',1:'orangered',2:'deepskyblue'}
# ,-1:'black'
cvr = colors.ColorConverter()
tmp = sorted(cols.keys())
cols_rgb = [cvr.to_rgb(cols[k]) for k in tmp]
intervals = np.asarray(tmp + [tmp[-1]+1]) - 0.5
cmap, norm = colors.from_levels_and_colors(intervals,cols_rgb)

legend_elements = [Patch(color='lightgrey', label=r'$n_{mut} < 3$'),
                    Patch(color='orangered', label=r'$P\nless P^{*}$'),
                    Patch(color='deepskyblue', label=r'$P < P^{*}$')]
# Patch(color='black', label="No data")

for taxon in pt.taxa:

    if taxon in pt.treatment_taxa_to_ignore:
        continue

    if taxon == 'J':
        continue

    if (taxon == 'C') or (taxon == 'P'):
        gene_name_fs = 10
    else:
        gene_name_fs = 14

    fig = plt.figure(figsize = (12, 20))

    ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=1)
    ax2 = plt.subplot2grid((1, 2), (0, 1), colspan=1)

    gene_names,gene_values = taxon_gene_dict[taxon]
    gene_names_zipped = list(zip(gene_names,gene_values))

    gene_names_zipped_1 = gene_names_zipped[:int(len(gene_names_zipped)/2)]
    gene_names_zipped_2 = gene_names_zipped[int(len(gene_names_zipped)/2):]


    gene_names_1 = [x[0] for x in gene_names_zipped_1]
    gene_values_1 = [x[1] for x in gene_names_zipped_1]
    gene_names_2 = [x[0] for x in gene_names_zipped_2]
    gene_values_2 = [x[1] for x in gene_names_zipped_2]

    data_1 = np.asarray(gene_values_1)
    data_2 = np.asarray(gene_values_2)

    ax1.xaxis.tick_top()
    ax1.set_xticks(np.arange(data_1.shape[1])+0.5, minor=False)
    ax1.set_yticks(np.arange(data_1.shape[0])+0.5, minor=False)
    ax1.set_yticklabels(gene_names_1, minor=False, fontsize=gene_name_fs)
    ax1.set_xticklabels(['1-Day', '10-Days', '100-Days'], minor=False, fontsize=14)

    ax2.xaxis.tick_top()
    ax2.set_xticks(np.arange(data_2.shape[1])+0.5, minor=False)
    ax2.set_yticks(np.arange(data_2.shape[0])+0.5, minor=False)
    ax2.set_yticklabels(gene_names_2, minor=False, fontsize=gene_name_fs)
    ax2.set_xticklabels(['1-Day', '10-Days', '100-Days'], minor=False, fontsize=14)


    ax1.legend(handles=legend_elements, bbox_to_anchor=(-0.65,1.12), loc="upper left", fontsize=16)
    #ax2.legend(handles=legend_elements, bbox_to_anchor=(-0.65,1.12), loc="upper left", fontsize=10)

    ax1.set_title(pt.latex_genus_bold_dict[taxon], weight='bold', fontsize=16)
    ax2.set_title(pt.latex_genus_bold_dict[taxon], weight='bold', fontsize=16)

    ax1.pcolor(data_1, cmap = cmap, norm = norm, edgecolors='k', linewidths=1.5)

    ax2.pcolor(data_2, cmap = cmap, norm = norm, edgecolors='k', linewidths=1.5)


    fig.subplots_adjust(wspace=.5, hspace=0.08)
    fig.savefig((pt.get_path() + '/figs/gene_convergence_plot_%s.pdf')%taxon, format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
