from __future__ import division
import os, sys, pickle, random
import numpy as np
from itertools import combinations

import  matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import ColorConverter

from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.stats as stats
import statsmodels.api as sm

import parse_file
import timecourse_utils
import mutation_spectrum_utils
import phylo_tools as pt

#import get_random_matrix

import phik

np.random.seed(123456789)
random.seed(123456789)

permutations_gene_content_divergence = 10000
#permutations_gene_content_divergence = 10

permutations_divergence = 10000

# permutations for anova
n_permutations = 100000
#n_permutations = 10
treatment_pairs = [['0','1'],['0','2'],['1','2']]


def calculate_IWR(array_1, array_2):

    # mutual information / joing entropy
    joint_entropy =  stats.entropy(array_1,array_2)



standardized_gene_overlap = {}
for taxon in pt.taxa:


    #if taxon == 'J':
    #    continue

    gene_dict = {}
    N_significant_genes_dict = {}

    gene_data = parse_file.parse_gene_list(taxon)
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    #locus_tag_to_gene_dict = {}
    #for gene_name_idx, gene_name in enumerate(gene_names):
    #    gene = genes[gene_name_idx]
    #    if gene == '':
    #        continue
    #    locus_tag_to_gene_dict[gene_name] = genes[gene_name_idx]


    if taxon == 'J':
        treatments_convergence = ['0', '2']

    else:
        treatments_convergence = ['0', '1', '2']


    for treatment in treatments_convergence:

        genes_significant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_%ss_%s.txt" % ('gene', treatment+taxon))
        genes_nonsignificant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_not_significant_%ss_%s.txt" % ('gene', treatment+taxon))

        if os.path.exists(genes_significant_file_path) == False:
            continue

        genes_significant_file = open(genes_significant_file_path, 'r')
        first_line_significant = genes_significant_file.readline()

        N_significant_genes = 0

        genes = []

        for line in genes_significant_file:
            line_split = line.strip().split(', ')
            gene_name = line_split[0]
            genes.append(gene_name)
            N_significant_genes += 1

        genes_significant_file.close()

        N_significant_genes_dict[treatment] = N_significant_genes
        gene_dict[treatment] = set(genes)


    for treatment_pair in combinations(treatments_convergence, 2):

        treatment_pair_set = set(treatment_pair)

        jaccard_treatment_pair = len(gene_dict[treatment_pair[0]] & gene_dict[treatment_pair[1]]) / len(gene_dict[treatment_pair[0]] | gene_dict[treatment_pair[1]])

        jaccard_null = []

        for i in range(permutations_gene_content_divergence):

            sample_1 = set(random.sample(gene_names, N_significant_genes_dict[treatment_pair[0]]))
            sample_2 = set(random.sample(gene_names, N_significant_genes_dict[treatment_pair[1]]))

            jaccard_treatment_pair_i = len(sample_1 & sample_2) / len(sample_1 | sample_2)

            jaccard_null.append(jaccard_treatment_pair_i)

        jaccard_null = np.asarray(jaccard_null)


        standardized_jaccard = (jaccard_treatment_pair-np.mean(jaccard_null)) / np.std(jaccard_null)

        if treatment_pair not in standardized_gene_overlap:
            standardized_gene_overlap[treatment_pair] = {}

        standardized_gene_overlap[treatment_pair][taxon] = {}
        standardized_gene_overlap[treatment_pair][taxon]['Z_jaccard'] = standardized_jaccard




        #standardized_gene_overlap[treatment_pair].append(standardized_jaccard)




# now do divergence


significant_multiplicity_dict = {}
significant_n_mut_dict = {}
gene_size_dict = {}
gene_mean_size_dict = {}
for taxon in pt.taxa:
    significant_multiplicity_dict[taxon] = {}
    significant_n_mut_dict[taxon] = {}
    gene_size_dict[taxon] = {}

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % ('0'+taxon)))
    Ltot = 0
    for gene_name in sorted(convergence_matrix.keys()):
        Lmin=0
        L = max([convergence_matrix[gene_name]['length'],Lmin])
        Ltot += L
    Lavg = Ltot*1.0/len(convergence_matrix.keys())

    gene_mean_size_dict[taxon] = Lavg

    for treatment_idx, treatment in enumerate(pt.treatments):

        significant_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significant_multiplicity_taxon_path) == False:
            continue
        significant_multiplicity_taxon = open(significant_multiplicity_taxon_path, "r")
        for i, line in enumerate( significant_multiplicity_taxon ):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            gene_size_dict[taxon][items[0]] = float(items[-5])
            if items[0] not in significant_multiplicity_dict[taxon]:
                significant_multiplicity_dict[taxon][items[0]] = {}

            if items[0] not in significant_n_mut_dict[taxon]:
                significant_n_mut_dict[taxon][items[0]] = {}

            significant_multiplicity_dict[taxon][items[0]][treatment] = float(items[-2])
            significant_n_mut_dict[taxon][items[0]][treatment] = float(items[-4])




sys.stderr.write("Runing permutational ANOVA for gene content....\n")


def gene_content_permutational_anova():


    #taxa_to_test = ['B','C','D','F','P','J']
    within_sum = 0
    all_divergenes_to_test = []
    all_mean_divergenes = []
    all_n_divergences = []
    for treatment_pair in standardized_gene_overlap.keys():

        if '1' in treatment_pair:
            taxa_to_test = ['B','C','D','F','P']
        else:
            taxa_to_test = pt.taxa

        divergences_treatment_pair = [standardized_gene_overlap[treatment_pair][taxon]['Z_jaccard'] for taxon in taxa_to_test]
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

        for taxon in pt.taxa:

            if taxon == 'J':

                #treatment_assignment = random.randrange(0,3)
                treatment_assignment = np.random.randint(0,3)
                J_div = standardized_gene_overlap[('0','2')][taxon]['Z_jaccard']

                if treatment_assignment == 0:
                    vs_1_10.append(J_div)
                elif treatment_assignment == 1:
                    vs_1_100.append(J_div)
                else:
                    vs_10_100.append(J_div)

            else:

                taxon_standardized_corr = []
                treatment_pairs_list = standardized_gene_overlap.keys()
                divergences_taxon = np.asarray([standardized_gene_overlap[l][taxon]['Z_jaccard'] for l in treatment_pairs_list])
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
    P_F = (len(F_permute_all[F_permute_all>F])+1) / (1+n_permutations)
    sys.stderr.write("Mean divergence: F = %.3f, P = %.4f \n" % (F, P_F))

    return F, P_F




F_gene, P_F_gene = gene_content_permutational_anova()


def calculate_divergence_correlations():

    sys.stdout.write("Starting divergence tests...\n")

    divergence_dict = {}

    for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

        treatment_pair_set = (treatment_pair[0], treatment_pair[1])

        divergence_dict[treatment_pair_set] = {}


        if '1' in treatment_pair:
            taxa = ['B','C','D','F','P']
        else:
            taxa = pt.taxa


        for taxon in taxa:

            #result = [(x[treatment_pair[0]],x[treatment_pair[1]]) for x in significant_multiplicity_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
            #result = [(x[treatment_pair[0]],x[treatment_pair[1]], x) for x in significant_n_mut_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
            result = [(dicts[treatment_pair[0]],dicts[treatment_pair[1]], keys) for keys, dicts in significant_n_mut_dict[taxon].items() if (treatment_pair[0] in dicts) and (treatment_pair[1] in dicts)]

            n_x = [int(x[0]) for x in result]
            n_y = [int(x[1]) for x in result]
            gene_names = [x[2] for x in result]

            gene_sizes_taxon_treatment_pair = [gene_size_dict[taxon][gene_i] for gene_i in gene_names]
            gene_sizes_taxon_treatment_pair = np.asarray(gene_sizes_taxon_treatment_pair)
            taxon_Lmean = gene_mean_size_dict[taxon]

            n_matrix = np.asarray([n_x, n_y])
            mult_matrix = n_matrix * (taxon_Lmean / gene_sizes_taxon_treatment_pair)
            rel_mult_matrix = mult_matrix/mult_matrix.sum(axis=1)[:,None]
            pearsons_corr = np.corrcoef(rel_mult_matrix[0,:], rel_mult_matrix[1,:])[1,0]
            pearsons_corr_squared = pearsons_corr**2

            pearsons_corr_squared_null = []
            for k in range(permutations_divergence):

                if (k % 2000 == 0) and (k>0):

                    sys.stdout.write("%d iterations\n" % (k))

                n_matrix_random = phik.simulation.sim_2d_data_patefield(n_matrix)
                mult_matrix_random = n_matrix_random * (taxon_Lmean / gene_sizes_taxon_treatment_pair)
                rel_mult_matrix_random = mult_matrix_random/mult_matrix_random.sum(axis=1)[:,None]
                pearsons_corr_random = np.corrcoef(rel_mult_matrix_random[0,:], rel_mult_matrix_random[1,:])[1,0]
                pearsons_corr_squared_random = pearsons_corr_random**2

                pearsons_corr_squared_null.append(pearsons_corr_squared_random)

            pearsons_corr_squared_null = np.asarray(pearsons_corr_squared_null)

            Z_corr = (pearsons_corr_squared - np.mean(pearsons_corr_squared_null)) / np.std(pearsons_corr_squared_null)

            P_corr = (len(pearsons_corr_squared_null[pearsons_corr_squared_null<pearsons_corr_squared])+1) / (permutations_divergence+1)

            divergence_dict[treatment_pair_set][taxon] = {}
            divergence_dict[treatment_pair_set][taxon]['pearsons_corr_squared'] = pearsons_corr_squared
            divergence_dict[treatment_pair_set][taxon]['P_value'] = P_corr
            divergence_dict[treatment_pair_set][taxon]['Z_corr'] = Z_corr

            sys.stdout.write("%d vs %d-day, %s: rho^2=%f, P=%f, Z=%f\n" % (10**int(treatment_pair[0]), 10**int(treatment_pair[1]), taxon, pearsons_corr_squared, P_corr, Z_corr))


    sys.stdout.write("Dumping pickle......\n")
    with open(pt.get_path()+'/data/divergence_pearsons.pickle', 'wb') as handle:
        pickle.dump(divergence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stdout.write("Done!\n")


#calculate_divergence_correlations()

with open(pt.get_path()+'/data/divergence_pearsons.pickle', 'rb') as handle:
    divergence_dict = pickle.load(handle)


sys.stderr.write("Runing permutational ANOVA for correlation....\n")


#taxa_to_test = ['B','C','D','F','P','J']
within_sum = 0
all_divergenes_to_test = []
all_mean_divergenes = []
all_n_divergences = []
for treatment_pair in divergence_dict.keys():

    if '1' in treatment_pair:
        taxa_to_test = ['B','C','D','F','P']
    else:
        taxa_to_test = pt.taxa

    divergences_treatment_pair = [divergence_dict[treatment_pair][taxon]['Z_corr'] for taxon in taxa_to_test]
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

    for taxon in pt.taxa:

        if taxon == 'J':

            #treatment_assignment = random.randrange(0,3)
            treatment_assignment = np.random.randint(0,3)
            J_div = divergence_dict[('0','2')][taxon]['Z_corr']

            if treatment_assignment == 0:
                vs_1_10.append(J_div)
            elif treatment_assignment == 1:
                vs_1_100.append(J_div)
            else:
                vs_10_100.append(J_div)

        else:

            taxon_standardized_corr = []
            treatment_pairs_list = divergence_dict.keys()
            divergences_taxon = np.asarray([divergence_dict[l][taxon]['Z_corr'] for l in treatment_pairs_list])
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



#from brokenaxes import brokenaxes

#ylim  = [28, 68]
#ylim2 = [-4, 4]
#ylimratio = (ylim[1]-ylim[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])
#ylim2ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])


gs = gridspec.GridSpec(nrows=2, ncols=1)

fig = plt.figure(figsize = (10, 13))
ax_divergence_gene = fig.add_subplot(gs[0, 0])
#fig, axes = plt.subplots(nrows=2, sharex=True)
#ax_divergence_gene_lower = axes[0]
#divider = make_axes_locatable(ax_divergence_gene_lower)
#ax_divergence_gene_upper = divider.new_vertical(size="100%", pad=0.1)
#fig.add_axes(ax_divergence_gene_upper)

#ax_divergence_gene_upper = fig.add_subplot(gs[0, 0])
#ax_divergence_gene_lower = fig.add_subplot(gs[1, 0])

#ylim_upper  = [30, 65]
#ylim_lower = [-1, 1]



#gs = gridspec.GridSpec(nrows=4, ncols=3)
#gs = gridspec.GridSpec(nrows=2, ncols=1)
#ax_divergence_gene = fig.add_subplot(gs[1, 0])

#ax_divergence_gene = brokenaxes(ylims=((-2, 3),(30, 65)), subplot_spec=gs[0, 0])


#fig, axes = plt.subplots(nrows=2, sharex=True)
#ax_divergence_gene = axes[0]
#ax_divergence_gene_lower = axes[0]
#divider = make_axes_locatable(ax_divergence_gene_lower)
#ax_divergence_gene_upper = divider.new_vertical(size="100%", pad=0.1)
#ax_divergence_gene_divider = divider.new_vertical(size="100%", pad=0.1)
#ax_divergence_gene = divider.new_vertical(size="100%", pad=0.1)

#fig.add_axes(ax_divergence_gene_upper)

#ylim_upper  = [30, 65]
#ylim_lower = [-1, 1]



#ax_divergence_gene = fig.add_subplot(gs[0, 0])
#ax_divergence_gene_break = fig.add_subplot(gs[1, 0], sharex=ax_divergence_gene)
ax_count=0
ax_count_divergence_gene=0

for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

    treatment_pair_set = (treatment_pair[0], treatment_pair[1])

    marker_style = dict(color='k', marker='o',
                markerfacecoloralt=pt.get_colors(treatment_pair[1]),
                markerfacecolor=pt.get_colors(treatment_pair[0]) )

    #standardized_gene_overlap_treatment_pair = standardized_gene_overlap[treatment_pair_set]

    if '1' in treatment_pair:
        taxa_to_test = ['B','C','D','F','P']
    else:
        taxa_to_test = pt.taxa

    standardized_gene_overlap_treatment_pair = [standardized_gene_overlap[treatment_pair_set][taxon]['Z_jaccard'] for taxon in taxa_to_test]

    divergence_gene_Z_mean = np.mean(standardized_gene_overlap_treatment_pair)
    divergence_gene_Z_se = np.std(standardized_gene_overlap_treatment_pair) / np.sqrt(len(standardized_gene_overlap_treatment_pair))

    ax_divergence_gene.errorbar(treatment_pair_idx, divergence_gene_Z_mean, yerr =divergence_gene_Z_se, \
            fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
            mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=2, ms=17)

    ax_divergence_gene.plot(treatment_pair_idx, np.mean(divergence_gene_Z_mean), markersize = 28,   \
        linewidth=2,  alpha=1, zorder=3, fillstyle='left', **marker_style)


    #ax_divergence_gene_upper.errorbar(treatment_pair_idx, divergence_gene_Z_mean, yerr =divergence_gene_Z_se, \
    #        fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
    #        mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=2, ms=17)


    #ax_divergence_gene_upper.plot(treatment_pair_idx, np.mean(divergence_gene_Z_mean), markersize = 25,   \
    #    linewidth=2,  alpha=1, zorder=3, fillstyle='left', **marker_style)


    #ax_divergence_gene_lower.errorbar(treatment_pair_idx, divergence_gene_Z_mean, yerr =divergence_gene_Z_se, \
    #        fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
    #        mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=2, ms=17)


    #ax_divergence_gene_lower.plot(treatment_pair_idx, np.mean(divergence_gene_Z_mean), markersize = 25,   \
    #    linewidth=2,  alpha=1, zorder=3, fillstyle='left', **marker_style)




#ax_divergence_gene.axhline( y=0, color='k', lw=3, linestyle=':', alpha = 1, zorder=1)

#ax_divergence.set_xticklabels([2,7.5,13], ['1-day vs. 10-days', '1-day vs. 100-days', '10-days vs. 100-days'], fontsize=13)
ax_divergence_gene.set_xticklabels( [], fontsize=12)

ax_divergence_gene.text(0.18, -0.04, '1-day vs. 10-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)
ax_divergence_gene.text(0.52, -0.04, '1-day vs. 100-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)
ax_divergence_gene.text(0.85, -0.04, '10-days vs. 100-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)

ax_divergence_gene.set_xlim([-0.5, 2.5])
ax_divergence_gene.set_ylim([30, 68])

#ax_divergence_gene.text(0.13, 0.945, 'Convergence', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)
#ax_divergence_gene.text(0.115, 0.83, 'Divergence', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)


ax_divergence_gene.tick_params(axis='x', labelsize=14, length = 0)

ax_divergence_gene.set_ylabel("Mean standardized Jaccard\nsimilarity of all taxa", fontsize = 16)

ax_divergence_gene.text(-0.05, 1.07, pt.sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)



ax_divergence_gene.set_title("Convergent/divergent evolution as the\nproportion of shared enriched genes ", fontsize=16, fontweight='bold')

ax_divergence_gene.arrow(0.06, 0.65, 0.0, 0.2, width=0.012,fc='k', ec='k', transform=ax_divergence_gene.transAxes)

ax_divergence_gene.text(0.06, 0.5, 'Increasing\nconvergence', fontsize=12, fontweight='bold', ha='center', va='center', rotation=90, transform=ax_divergence_gene.transAxes)


ax_divergence_gene.text(0.865, 0.17, r'$F=%s$' % "{0:.3g}".format(F_gene), fontsize=14, ha='center', va='center', transform=ax_divergence_gene.transAxes)
ax_divergence_gene.text(0.865, 0.1, r'$P=%s$' % "{0:.3g}".format(P_F_gene), fontsize=14, ha='center', va='center', transform=ax_divergence_gene.transAxes)




#ax_divergence_gene.plot()


#ax_divergence_gene_upper.set_ylim(ylim)
#ax_divergence_gene_lower.set_ylim(ylim2)

#ax_divergence_gene_upper.spines['bottom'].set_visible(False)
#ax_divergence_gene_lower.spines['top'].set_visible(False)
#ax_divergence_gene_upper.xaxis.tick_top()
#ax_divergence_gene_upper.tick_params(labeltop='off')
#ax_divergence_gene_lower.xaxis.tick_bottom()



#ax_divergence_gene_upper.set_ylim(ylim_upper)
#ax_divergence_gene_lower.set_ylim(ylim_lower)

#kwargs = dict(color='k', clip_on=False)
#xlim = ax_divergence_gene_upper.get_xlim()
#dx = .02*(xlim[1]-xlim[0])
#dy = .01*(ylim[1]-ylim[0])/ylimratio
#ax_divergence_gene_upper.plot((xlim[0]-dx,xlim[0]+dx), (ylim[0]-dy,ylim[0]+dy), **kwargs)
#ax_divergence_gene_upper.plot((xlim[1]-dx,xlim[1]+dx), (ylim[0]-dy,ylim[0]+dy), **kwargs)
#dy = .01*(ylim2[1]-ylim2[0])/ylim2ratio
#ax_divergence_gene_lower.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
#ax_divergence_gene_lower.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
#ax_divergence_gene_upper.set_xlim(xlim)
#ax_divergence_gene_lower.set_xlim(xlim)


#ax_divergence_gene_upper.spines['bottom'].set_visible(False)
#ax_divergence_gene_lower.spines['top'].set_visible(False)






#axis_break1 = 2
#axis_break2 = 30
#x_min = -0.75
#x_max = len(list(enumerate(treatment_pairs))) +0.5
#l = 0.2  # "break" line length
#kwargs = dict(color="k", clip_on=False, linewidth=1)
#ax_divergence_gene_upper.plot((x_min - l, x_min + l), (axis_break2, axis_break2), **kwargs)# top-left
#ax_divergence_gene_upper.plot((x_max - l, x_max + l), (axis_break2, axis_break2), **kwargs)# top-right
#ax_divergence_gene_lower.plot((x_min - l, x_min + l), (axis_break1, axis_break1), **kwargs)# bottom-left
#ax_divergence_gene_lower.plot((x_max - l, x_max + l), (axis_break1, axis_break1), **kwargs)# bottom-right


#ylim_ratio_upper = (ylim_upper[1]-ylim_upper[0])/(ylim_lower[1]-ylim_lower[0]+ylim_upper[1]-ylim_upper[0])
#ylim_ratio_lower = (ylim_lower[1]-ylim_lower[0])/(ylim_lower[1]-ylim_lower[0]+ylim_upper[1]-ylim_upper[0])

#kwargs = dict(color='k', clip_on=False)
#xlim = ax_divergence_gene_upper.get_xlim()
#dx = .02*(xlim[1]-xlim[0])
#dy = .01*(ylim_upper[1]-ylim_upper[0])/ylim_ratio_upper
#ax_divergence_gene_upper.plot((xlim[0]-dx,xlim[0]+dx), (ylim_upper[0]-dy,ylim_upper[0]+dy), **kwargs)
#ax_divergence_gene_upper.plot((xlim[1]-dx,xlim[1]+dx), (ylim_upper[0]-dy,ylim_upper[0]+dy), **kwargs)
#dy = .01*(ylim_lower[1]-ylim_lower[0])/ylim_ratio_lower
#ax_divergence_gene_lower.plot((xlim[0]-dx,xlim[0]+dx), (ylim_lower[1]-dy,ylim_lower[1]+dy), **kwargs)
#ax_divergence_gene_lower.plot((xlim[1]-dx,xlim[1]+dx), (ylim_lower[1]-dy,ylim_lower[1]+dy), **kwargs)

#ax_divergence_gene_upper.set_xlim(xlim)
#ax_divergence_gene_lower.set_xlim(xlim)

# From https://matplotlib.org/examples/pylab_examples/broken_axis.html
#d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them


#d = .25  # proportion of vertical to horizontal extent of the slanted line
#kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
#              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
#ax_divergence_gene_upper.plot([0, 1], [0, 0], transform=ax_divergence_gene_upper.transAxes, **kwargs)
#ax_divergence_gene_lower.plot([0, 1], [1, 1], transform=ax_divergence_gene_lower.transAxes, **kwargs)


#
#break_y_axis = 0.5
#kwargs = dict(transform=ax_divergence_gene_upper.transAxes, color='k', clip_on=False)
#ax_divergence_gene_upper.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
#ax_divergence_gene_upper.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

#kwargs.update(transform=ax_divergence_gene_lower.transAxes)  # switch to the bottom axes
#ax_divergence_gene_lower.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
#ax_divergence_gene_lower.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal



#ax_divergence_gene.spines['bottom'].set_visible(False)
#ax_divergence_gene.xaxis.tick_top()
#ax_divergence_gene.tick_params(labeltop='off')  # don't put tick labels at the top

#ax_divergence_gene_break.spines['top'].set_visible(False)
#ax_divergence_gene_break.xaxis.tick_bottom()


#ax_divergence_gene.set_xticklabels( [], fontsize=12)

#ax_divergence_gene.text(0.15, -0.04, '1-day vs. 10-days', fontsize=11, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)
#ax_divergence_gene.text(0.5, -0.04, '1-day vs. 100-days', fontsize=11, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)
#ax_divergence_gene.text(0.84, -0.04, '10-days vs. 100-days', fontsize=11, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)

#ax_divergence_gene.set_xlim([-0.5, 2.5])
#ax_divergence_gene.set_ylim([-12, 1.5])

#ax_divergence_gene_lower.text(0.13, 0.945, 'Convergence', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene_lower.transAxes)
#ax_divergence_gene_lower.text(0.115, 0.83, 'Divergence', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene_lower.transAxes)

#ax_divergence_gene_lower.axhline( y=0, color='k', lw=3, linestyle=':', alpha = 1, zorder=1)

#ax_divergence_gene.tick_params(axis='x', labelsize=14, length = 0)

#ax_divergence_gene.set_ylabel("Mean standardized Jaccard index of enriched genes, "+ r'$\bar{Z}_{\rho^{2}}$' , fontsize = 15)

#ax_divergence_gene.text(-0.05, 1.07, pt.sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_divergence_gene.transAxes)





#ax_divergence = fig.add_subplot(gs[2:4, 0:3])
ax_divergence = fig.add_subplot(gs[1, 0])
#ax_divergence = axes[1]
ax_divergence.axhline( y=0, color='k', lw=3, linestyle=':', alpha = 1, zorder=1)
ax_count_divergence=0
ax_count+=1
for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

    treatment_pair_set = (treatment_pair[0], treatment_pair[1])



    ax_count_divergence_treatment_pair = []
    #treatment_pair_slopes = []

    if '1' in treatment_pair:
        taxa = ['B','C','D','F','P']
    else:
        taxa = pt.taxa

    divergence_Z_pearsons_pair = []


    for taxon in taxa:

        standardized_correlation = divergence_dict[treatment_pair_set][taxon]['Z_corr']

        divergence_Z_pearsons_pair.append(standardized_correlation)


        ax_count_divergence+=1

    marker_style = dict(color='k', marker='o',
                markerfacecoloralt=pt.get_colors(treatment_pair[1]),
                markerfacecolor=pt.get_colors(treatment_pair[0]) )

    divergence_Z_mean = np.mean(divergence_Z_pearsons_pair)
    divergence_Z_se = np.std(divergence_Z_pearsons_pair) / np.sqrt(len(divergence_Z_pearsons_pair))

    ax_divergence.errorbar(treatment_pair_idx, divergence_Z_mean, yerr =divergence_Z_se, \
            fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
            mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=2, ms=17)


    ax_divergence.plot(treatment_pair_idx, np.mean(divergence_Z_pearsons_pair), markersize = 28,   \
        linewidth=2,  alpha=1, zorder=3, fillstyle='left', **marker_style)



    #if treatment_pair_idx < 2:

    #    ax_divergence.axvline( x=ax_count_divergence-0.5, color='k', lw=2, linestyle='-', alpha = 1, zorder=2)







ax_divergence.text(0.865, 0.17, r'$F=%s$' % "{0:.3g}".format(F), fontsize=14, ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.865, 0.1, r'$P=%s$' % "{0:.3g}".format(P_F), fontsize=14, ha='center', va='center', transform=ax_divergence.transAxes)





#ax_divergence.set_xticklabels([2,7.5,13], ['1-day vs. 10-days', '1-day vs. 100-days', '10-days vs. 100-days'], fontsize=13)
ax_divergence.set_xticklabels( [], fontsize=12)

ax_divergence.text(0.18, -0.04, '1-day vs. 10-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.52, -0.04, '1-day vs. 100-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.85, -0.04, '10-days vs. 100-days', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)



ax_divergence.set_xlim([-0.5, 2.5])
ax_divergence.set_ylim([-12, 1.5])

ax_divergence.text(0.13, 0.925, 'Convergence', fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.115, 0.85, 'Divergence', fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)


ax_divergence.tick_params(axis='x', labelsize=14, length = 0)

ax_divergence.set_ylabel("Mean standardized\ncorrelation of all taxa, "+ r'$\bar{Z}_{\rho^{2}}$' , fontsize = 16)

ax_divergence.text(-0.05, 1.07, pt.sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)


ax_divergence.set_title("Convergent/divergent evolution as the\ncorrelation in mutation counts across enriched genes", fontsize=16, fontweight='bold')

#ax_divergence_gene.set_title("Convergent/divergent evolution as gene identity", fontsize=15, fontweight='bold')

fig.subplots_adjust(hspace=0.25,wspace=0.2) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/per_taxon _divergence.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



sys.stderr.write("Done with figure!\n")
