from __future__ import division
import os, sys, itertools
import numpy as np

import phylo_tools as pt
import parse_file
import timecourse_utils
from itertools import combinations

import  matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

import scipy.stats as stats

import statsmodels.stats.multitest as multitest

np.random.seed(123456789)

sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']


set_time_dict = {'B':500, 'C':500, 'D':500, 'P':600, 'F':500, 'J':500}

mutation_trajectories = {}

#taxa = pt.taxa
taxa = ['B','C','D','P', 'F', 'J']
# loop through taxa and get M(700) for all reps in each treatment
for taxon in taxa:

    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment in treatments:
        for replicate in pt.replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            mutation_trajectories[population] = {}

            #sys.stderr.write("Processing %s...\t" % population)

            times, Ms, fixed_Ms = parse_file.get_mutation_fixation_trajectories(population)


            #time_idx = np.where(times == set_time)
            for time in times:
                time_idx = np.where(times == time)

                mutation_trajectories[population][time] = (times[time_idx][0],np.log10(Ms[time_idx][0] /  time))


            #sys.stderr.write("analyzed %d mutations in %s!\n" % (len(Ms) ,population))


#fig = plt.figure(figsize = (16, 4)) #
#fig.subplots_adjust(bottom= 0.15)

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(nrows=6, ncols=4)


legend_elements = [Line2D([0], [0], c='k', lw=2, linestyle='--', label='Estimated ' + r'$\beta_{1}$'),
                   Line2D([0], [0], c='darkgrey', linestyle=':', lw=2, label='Null ' + r'$\beta_{1}$')]


#gs = gridspec.GridSpec(nrows=3, ncols=4)
#gs = gridspec.GridSpec(nrows=6, ncols=5)

x_axis_labels = ['Tra nsfer time (days)']
#slope_null = [-1,1]
slope_null = -1
sub_plot_counts = 0

p_values = []
taxon_axis_dict = {}

text_position_dict = {'B':(0.6,0.9), 'C':(0.6,0.9),  'D':(0.6 ,0.9), 'F':(0.4, 0.9), 'P':(0.4, 0.9)}

for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):

    for taxon_idx, taxon in enumerate(taxon_list):

        if taxon == 'J':
            treatments = ['0','2']
        else:
            treatments = pt.treatments

        set_time = set_time_dict[taxon]

        #ax = fig.add_subplot(gs[taxon_idx, taxon_list_idx])
        ax = fig.add_subplot(gs[taxon_idx*2:(taxon_idx*2)+2  , taxon_list_idx])


        taxon_axis_dict[taxon] = ax

        #if time_measure_idx == 0:
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )

        ax.set_xlabel('Transfer time (days)', fontsize=10)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.xaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_tick_params(labelsize=8)

        #ax.set_ylabel('Mutations per-day,' + str(set_time) + ', ' + r'$M({{{}}})$'.format(set_time), fontsize=10)

        ax.set_ylabel('Mutations per-day, $M(t)/t$', fontsize=11)

        ax.text(-0.1, 1.07, sub_plot_labels[sub_plot_counts], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

        if sub_plot_counts == 0:
            ax.legend(handles=legend_elements, loc='lower left',  prop={'size': 7})
        # residuals plots

        sub_plot_counts += 1

        times_all_list = []
        mutations_all_list = []

        for treatment in treatments:

            #print(mutation_trajectories.items())
            mutations_list = []
            for key, value in mutation_trajectories.items():
                if (treatment+taxon in key) and (set_time in value):

                    mutations_list.append(value[set_time][1])

            if len(mutations_list) == 0:
                continue

            mutations_list = np.asarray(mutations_list) #/ set_time_dict[taxon]
            #mutations_list = np.asarray([value[set_time][1] for key, value in mutation_trajectories.items() if (treatment+taxon in key) and (set_time in value.values())])
            times_list = np.repeat( int(treatment), len(mutations_list))

            ax.scatter((10**times_list) + np.random.randn(len(times_list))*0.1 , 10**mutations_list, s= 140, linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=0.8, zorder=3)
            times_all_list.extend(times_list)
            mutations_all_list.extend(mutations_list)


        ax.set_ylim([(10** min(mutations_all_list))*0.5, (10** max(mutations_all_list))*2  ])
        ax.set_xlim([(10** min(times_all_list))*0.5, (10** max(times_all_list))*2  ])

        mutations_all_list = np.asarray(mutations_all_list ) #/ set_time_dict[taxon]


        if taxon == 'J':
            treatment_1 = [x[1] for x in list(zip(times_all_list, mutations_all_list)) if x[0] ==0]
            treatment_2 = [x[1] for x in list(zip(times_all_list, mutations_all_list)) if x[0] ==2]

            t, p_value = stats.ttest_ind(treatment_1, treatment_2)
            ax.text(0.5, 0.9, r'$t=$' + str(round(t, 2)), fontsize=9, transform=ax.transAxes)

            #if p_value < 0.05:
            #    ax.text(0.5, 0.8, r'$P < 0.05$', fontsize=9, transform=ax.transAxes)
            #else:
            #    ax.text(0.5, 0.8, r'$P \nless 0.05$', fontsize=9, transform=ax.transAxes)
            ax.text(0.5, 0.8, r'$P = $'+ str(round(p_value, 3)), fontsize=9, transform=ax.transAxes)


        else:

            slope, intercept, r_value, p_value, std_err = stats.linregress(times_all_list, mutations_all_list)
            x_log10_fit_range =  np.linspace(min(times_all_list) * 0.5, max(times_all_list) * 1.5, 10000)


            y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)
            ax.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', zorder=2)

            y_null_range = 10 ** (slope_null * x_log10_fit_range + intercept)
            ax.plot(10**x_log10_fit_range, y_null_range, c='darkgrey', linestyle=':', lw=3, zorder=1)

            # hypothetical slope of -1
            ratio = (slope - slope_null) / std_err
            pval = stats.t.sf(np.abs(ratio), len(mutations_all_list)-2)#*2
            #SS_residual =

            p_values.append(pval)
            # one tailed
            sys.stderr.write("Species %s slope t-test = %f, p = %f\n" % (taxon , round(ratio, 3), round(pval, 3)))

            taxon_x = text_position_dict[taxon][0]
            taxon_y = text_position_dict[taxon][1]

            ax.text(taxon_x, taxon_y, r'$\beta_{1}=$' + str(round(slope, 2)), fontsize=9, transform=ax.transAxes)


            #ax.text(0.05, 0.13, r'$r^{2}=$' + str(round(r_value**2, 2)), fontsize=8, transform=ax.transAxes)



reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')


for taxon_idx, taxon in enumerate(['B','C','D', 'F','P']):

    if taxon == 'J':
        continue

    p_value = pvals_corrected[taxon_idx]

    taxon_x = text_position_dict[taxon][0]
    taxon_y = text_position_dict[taxon][1] - 0.1
    #if p_value < 0.05:
    #    taxon_axis_dict[taxon].text(taxon_x, taxon_y, r'$P_{BH} < 0.05$', fontsize=9, transform=taxon_axis_dict[taxon].transAxes)
    #else:
    #    taxon_axis_dict[taxon].text(taxon_x, taxon_y, r'$P_{BH} \nless 0.05$', fontsize=9, transform=taxon_axis_dict[taxon].transAxes)

    #taxon_axis_dict[taxon].text(taxon_x, taxon_y, r'$P_{BH} < $' + str(round(p_value, 3)), fontsize=9, transform=taxon_axis_dict[taxon].transAxes)
    taxon_axis_dict[taxon].text(taxon_x, taxon_y, r'$P_{{BH}} < 10^{{{}}}$'.format(str( int(np.log10(p_value)) )), fontsize=9, transform=taxon_axis_dict[taxon].transAxes)





# then plot fixed mutations
fixed_mutations_per_day = {}
for taxon in taxa:
    if taxon == 'J':
        treatments = ['0','2']
    else:
        treatments = pt.treatments

    for treatment in treatments:
        for replicate in pt.replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            times, Ms, fixed_Ms = parse_file.get_mutation_fixation_trajectories(population)

            if isinstance(fixed_Ms,float) == True:
                fixed_Ms = np.asarray([0]* len(times))

            #fixed_mutations_per_day[population] = (fixed_Ms[-1]+(int(times[-1])/1000))/int(times[-1])
            fixed_mutations_per_day[population] = fixed_Ms[-1]/times[-1]


ax_pca = fig.add_subplot(gs[0:3, 2:4])
#gs[taxon_idx*2:(taxon_idx*2)+2  , taxon_list_idx]

count = 0
for treatment in pt.treatments:

    for taxon in ['B','C','D','F','J','P']:

        if (taxon == 'J') and (treatment == '1'):
            continue

        fixed_mutations_i = [y for x,y in fixed_mutations_per_day.items() if treatment+taxon in x]
        #fixed_mean = np.mean(np.log10(fixed_mutations_i))
        #fixed_std = np.std(np.log10(fixed_mutations_i))

        fixed_mean = np.mean(fixed_mutations_i)
        fixed_std = np.std(fixed_mutations_i)

        #ax_pca.errorbar(fixed_mean, count, xerr=2*fixed_std, linestyle='-', marker='o', lw = 3)


        #ax_pca.errorbar(fixed_mean, count, xerr=(2*fixed_std), \
        #        fmt = 'o', alpha = 1, barsabove = True, marker = 's', \
        #        mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=1, ms=17)


        #ax_pca.scatter(fixed_mutations_i, list(itertools.repeat(count, len(fixed_mutations_i))), s = 250, \
        #    linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), \
        #    edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), \
        #    alpha=0.8, zorder=2)

        ax_pca.scatter(list(itertools.repeat(count, len(fixed_mutations_i))), fixed_mutations_i, s = 180, \
            linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), \
            edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), \
            alpha=0.6, zorder=2)

        count+=1


ax_pca.set_ylabel('Fixed mutations, per-day', fontsize=14)
ax_pca.text(-0.1, 1.01, sub_plot_labels[sub_plot_counts], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_pca.transAxes)

#ax_pca.set_xscale('symlog')

ax_pca.set_xticks([])
# for minor ticks
ax_pca.set_xticks([], minor=True)
ax_pca.set_xticklabels([])


ax_pca.set_xticks([3, 8, 13])
ax_pca.set_xticklabels(['1-day', '10-days', '100-days'],fontweight='bold' )
ax_pca.tick_params(axis='x', labelsize=12, length = 0)

#ax_pca.text(-0.05, 0.2, '1-day', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)
#ax_pca.text(-0.05, 0.5, '10-days', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)
#ax_pca.text(-0.05, 0.8, '100-days', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)

sub_plot_counts += 1




# plot fmax curve


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


ax_fmax = fig.add_subplot(gs[3:6, 2:4])

for treatment in pt.treatments:

    for taxon, f_max_array in fmax_dict[treatment].items():

        f_max_array_sort = np.sort(f_max_array)
        cdf = 1-  np.arange(len(f_max_array_sort))/float(len(f_max_array_sort))
        #num_bins = 40
        #counts, bin_edges = np.histogram(f_max_array_sort, bins=num_bins, normed=True)
        #cdf = np.cumsum(counts)
        #pylab.plot(bin_edges[1:], cdf)

        ax_fmax.plot(f_max_array_sort, cdf, c =pt.get_colors(treatment), alpha=0.8)
        #marker=pt.plot_species_marker(taxon), markersize=1)



ax_fmax.set_xlim([ 0.09, 1.03 ])
ax_fmax.set_ylim([ 0.0008, 1.03 ])

ax_fmax.set_xscale('log', base=10)
ax_fmax.set_yscale('log', base=10)

ax_fmax.set_xlabel('Maximum observed allele frequency, ' + r'$f_{max}$', fontsize=12)
ax_fmax.set_ylabel('Fraction of mutations ' + r'$\geq f_{max}$', fontsize=13)
ax_fmax.text(-0.1, 1.01, sub_plot_labels[sub_plot_counts], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_fmax.transAxes)


ins_ks = inset_axes(ax_fmax, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.12,0.07,0.4,0.38), bbox_transform=ax_fmax.transAxes)
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


    mean_D_dict[treatment_pair] = {}
    mean_D = np.mean(y)
    se_D = np.std(y) / np.sqrt(len(y))
    mean_D_dict[treatment_pair]['mean_D'] = mean_D
    mean_D_dict[treatment_pair]['se_D'] = se_D

    marker_style = dict(color='k',
                        markerfacecoloralt=pt.get_colors(treatment_pair[0]),
                        markerfacecolor=pt.get_colors(treatment_pair[1]))

    if treatment_pair_idx == 0:
        count_real = count + 0.3
    elif treatment_pair_idx == 1:
        count_real = count
    else:
        count_real = count - 0.3

    ins_ks.errorbar(count_real, mean_D, yerr = se_D , \
            fmt = 'o', alpha = 1, barsabove = True, marker = 'o', \
            mfc = 'white', mec = 'white', lw=1.5, c = 'k', zorder=1, ms=17)

    ins_ks.plot(count_real, mean_D, markersize = 9, marker = 'o',  \
        linewidth=0.2,  alpha=1, fillstyle='left', zorder=2, **marker_style)

    count+=1


ins_ks.tick_params(labelsize=5)
ins_ks.tick_params(axis='both', which='major', pad=1)

ins_ks.set_xticks([0.3, 1, 1.7])
ins_ks.set_xticklabels(['1 vs.\n10-days', '1 vs.\n100-days', '10 vs.\n100-days'],fontweight='bold' )
ins_ks.tick_params(axis='x', labelsize=6.5, length = 0)

#temp = ins_ks.yaxis.get_ticklabels()
#temp = list(set(temp) - set(temp[::2]))
#for label in temp:
#    label.set_visible(False)


ins_ks.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


ins_ks.set_xlim([-0.04,1.92])
ins_ks.set_ylim([0.05,0.25])


# permutational anova




sys.stderr.write("Runing permutational ANOVA for f_max distance....\n")

n_permutations=10000

#taxa_to_test = ['B','C','D','F','P','J']
within_sum = 0
all_divergenes_to_test = []
all_mean_divergenes = []
all_n_divergences = []
for treatment_pair in ks_dict.keys():

    if '1' in treatment_pair:
        taxa_to_test = ['B','C','D','F','P']
    else:
        taxa_to_test = pt.taxa

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

    for taxon in pt.taxa:

        if taxon == 'J':

            #treatment_assignment = random.randrange(0,3)
            treatment_assignment = np.random.randint(0,3)
            J_div = ks_dict[('0','2')][taxon]['D']

            if treatment_assignment == 0:
                vs_1_10.append(J_div)
            elif treatment_assignment == 1:
                vs_1_100.append(J_div)
            else:
                vs_10_100.append(J_div)

        else:

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


ins_ks.set_ylim([0.05,0.25])
ins_ks.text(0.3, 0.23, r'$F=%s$' % "{0:.3g}".format(F), fontsize=9, ha='center', va='center', transform=ins_ks.transAxes)
ins_ks.text(0.35, 0.1, r'$P=%s0$' % "{0:.4g}".format(P_F), fontsize=9, ha='center', va='center', transform=ins_ks.transAxes)


sys.stderr.write("Saving figure....\n")

fig.subplots_adjust(hspace=1.8, wspace=0.5) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/mutation_regression.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
