from __future__ import division
import os, sys, itertools
import numpy as np

import phylo_tools as pt
import parse_file
import timecourse_utils

import  matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import scipy.stats as stats

import statsmodels.stats.multitest as multitest


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

fig = plt.figure(figsize = (12, 9))
gs = gridspec.GridSpec(nrows=3, ncols=4)

x_axis_labels = ['Transfer time (days)']
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

        ax = fig.add_subplot(gs[taxon_idx, taxon_list_idx])
        taxon_axis_dict[taxon] = ax

        #if time_measure_idx == 0:
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )

        ax.set_xlabel('Transfer time (days)', fontsize=11)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.xaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_tick_params(labelsize=8)

        #ax.set_ylabel('Mutations per-day,' + str(set_time) + ', ' + r'$M({{{}}})$'.format(set_time), fontsize=10)

        ax.set_ylabel('Mutations per-day, $M(t)/t$', fontsize=11)

        ax.text(-0.1, 1.07, sub_plot_labels[sub_plot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

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
            ax.text(0.5, 0.9, r'$t=$' + str(round(t, 2)), fontsize=8, transform=ax.transAxes)

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


        ax_pca.scatter(fixed_mutations_i, list(itertools.repeat(count, len(fixed_mutations_i))), s = 250, \
            linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), \
            edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), \
            alpha=0.8, zorder=2)

        count+=1


ax_pca.set_xlabel('Fixed mutations, per-day', fontsize=14)
ax_pca.text(-0.1, 1.01, sub_plot_labels[sub_plot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_pca.transAxes)

#ax_pca.set_xscale('symlog')

ax_pca.set_yticks([])
# for minor ticks
ax_pca.set_yticks([], minor=True)
ax_pca.set_yticklabels([])

ax_pca.text(-0.05, 0.2, '1-day', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)
ax_pca.text(-0.05, 0.5, '10-days', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)
ax_pca.text(-0.05, 0.8, '100-days', ha='center', va='center', rotation='vertical', fontsize=16, transform=ax_pca.transAxes)


fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/mutation_regression.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
