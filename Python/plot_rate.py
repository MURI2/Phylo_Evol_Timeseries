from __future__ import division
import os, sys
import numpy as np

import phylo_tools as pt
import parse_file
import timecourse_utils

import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D


import scipy.stats as stats




def plot_mutation_trajectory_taxon(taxon):

    if taxon == 'J':
        treatments = ['0','2']
        sub_plot_labels = ['a','b', 'c','d', 'e','f', 'g','h', 'i','j']
        sub_plot_count_step = 2
        dim = (6, 15)
    else:
        treatments = pt.treatments
        sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i', 'j','k','l', 'm','n','o']
        sub_plot_count_step = 3
        dim=(10, 15)


    sys.stderr.write("Loading mutation data...\n")

    mutation_trajectories = {}
    fixed_mutation_trajectories = {}
    delta_mutation_trajectories = {}
    #transit_times = {}
    median_trajectories = {}
    n_muts_trajectories = {}

    for treatment in treatments:
        for replicate in pt.replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            sys.stderr.write("Processing %s...\t" % population)

            times, Ms, fixed_Ms = parse_file.get_mutation_fixation_trajectories(population)

            times_, medians_log10, num_muts = parse_file.get_mutation_fixation_trajectories_median_freq_and_mut_number(population)

            if isinstance(fixed_Ms,float) == True:
                fixed_Ms = np.asarray([0]* len(times))

            fixed_mutation_trajectories[population] = (times, fixed_Ms)
            mutation_trajectories[population] = (times,np.log10(Ms))
            delta_mutation_trajectories[population] = (times[1:], np.log10(Ms[1:]/Ms[:-1] ))

            median_trajectories[population] = (times_, medians_log10)
            n_muts_trajectories[population] = (times_, num_muts)

            sys.stderr.write("analyzed %d mutations!\n" % len(Ms))

    fig = plt.figure(figsize = dim)

    column_count = 0

    for treatment in treatments:

        ax_t_vs_M = plt.subplot2grid((5, len(treatments)), (0, column_count), colspan=1)
        ax_t_vs_delta_M = plt.subplot2grid((5, len(treatments)), (1, column_count), colspan=1)
        ax_t_vs_F = plt.subplot2grid((5, len(treatments)), (2, column_count), colspan=1)

        ax_t_vs_median_freq = plt.subplot2grid((5, len(treatments)), (3, column_count), colspan=1)
        ax_t_vs_number_muts = plt.subplot2grid((5, len(treatments)), (4, column_count), colspan=1)

        ax_t_vs_M.text(-0.1, 1.07, sub_plot_labels[column_count], fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_t_vs_M.transAxes)
        ax_t_vs_delta_M.text(-0.1, 1.07, sub_plot_labels[column_count+sub_plot_count_step], fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_t_vs_delta_M.transAxes)
        ax_t_vs_F.text(-0.1, 1.07, sub_plot_labels[column_count+sub_plot_count_step*2], fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_t_vs_F.transAxes)
        ax_t_vs_median_freq.text(-0.1, 1.07, sub_plot_labels[column_count+sub_plot_count_step*3], fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_t_vs_median_freq.transAxes)
        ax_t_vs_number_muts.text(-0.1, 1.07, sub_plot_labels[column_count+sub_plot_count_step*4], fontsize=14, fontweight='bold', ha='center', va='center', transform=ax_t_vs_number_muts.transAxes)


        treatment_taxon_populations = []


        all_medians = []
        all_numbers = []

        for replicate in pt.replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            Mts,Ms = mutation_trajectories[population]
            fixed_Mts, fixed_Ms = fixed_mutation_trajectories[population]
            deltaMts, deltaMs = delta_mutation_trajectories[population]


            median_trajectories_ts, median_trajectories_ = median_trajectories[population]
            n_muts_trajectories_ts, n_muts_trajectories_ = n_muts_trajectories[population]


            ax_t_vs_M.plot(Mts, 10**Ms, 'o-',color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            ax_t_vs_M.set_yscale('log', base=10)
            ax_t_vs_M.tick_params(axis='x', labelsize=8)

            # back transform to format plot axes
            ax_t_vs_delta_M.plot(deltaMts, 10**deltaMs, color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon))
            ax_t_vs_delta_M.set_yscale('log', base=10)

            ax_t_vs_F.plot(fixed_Mts, fixed_Ms, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            #ax_M_vs_F.set_xlabel('Days, ' + r'$t$', fontsize = 12)




            ax_t_vs_median_freq.plot(median_trajectories_ts, 10**median_trajectories_, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            ax_t_vs_median_freq.set_yscale('log', base=10)
            #ax_t_vs_median_freq.tick_params(axis='y', labelsize=6)

            ax_t_vs_median_freq.yaxis.set_tick_params(labelsize=8)

            all_medians.extend(median_trajectories_.tolist())



            ax_t_vs_number_muts.plot(n_muts_trajectories_ts, n_muts_trajectories_, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            ax_t_vs_number_muts.set_yscale('log', base=10)

            ax_t_vs_number_muts.tick_params(axis='y', labelsize=8)

            all_numbers.extend(n_muts_trajectories_.tolist())



            treatment_taxon_populations.append(population)


        print(10**(min(all_medians)*0.8), 10**( max(all_medians)*1.2))

        ax_t_vs_median_freq.set_ylim([10**(min(all_medians))*0.8, 10**( max(all_medians))*1.2])


        ax_t_vs_number_muts.set_ylim([min(all_numbers)*0.8, max(all_numbers)*1.2])


        avg_Mts, avg_Ms = timecourse_utils.average_trajectories([mutation_trajectories[population] for population in treatment_taxon_populations])

        avg_deltaMts, avg_deltaMs = timecourse_utils.average_trajectories([delta_mutation_trajectories[population] for population in treatment_taxon_populations])

        ax_t_vs_delta_M.axhline(y=1, c='grey', linestyle=':', lw=3, zorder=1)
        ax_t_vs_M.plot(avg_Mts, 10**avg_Ms, '--',color='k', marker=" ", alpha=1, linewidth=4, zorder=2)
        ax_t_vs_delta_M.plot(avg_deltaMts, 10**avg_deltaMs, '--',color='k', marker=" ", alpha=1, linewidth=4, zorder=2)

        # keep them on the same y axes
        if taxon == 'C':
            ax_t_vs_delta_M.set_ylim([0.2,42])
        elif taxon == 'D':
            ax_t_vs_delta_M.set_ylim([0.2,20])

        if (column_count==0):
            legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}(t)$')]
            ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=8)

        ax_t_vs_M.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

        #if treatment == '2':
        #    ax_M_vs_F.yaxis.set_major_locator(MaxNLocator(integer=True))

        if column_count == 0:

            ax_t_vs_M.set_ylabel('Mutations, ' + r'$M(t)$', fontsize = 15)
            ax_t_vs_F.set_ylabel('Fixed mutations', fontsize = 15)
            ax_t_vs_delta_M.set_ylabel('Change in mutations,\n' + r'$M(t)/M(t-1)$', fontsize = 15)

            ax_t_vs_median_freq.set_ylabel('Median mutation freq.\nat time $t$', fontsize = 15)
            ax_t_vs_number_muts.set_ylabel('Number of mutations\nat time $t$', fontsize = 15)


        column_count += 1

    fig.text(0.53, 0.05, 'Days, ' + r'$t$', ha='center', fontsize=28)
    fig.suptitle(pt.latex_genus_dict[taxon], fontsize=30)
    fig_name = pt.get_path() + '/figs/rate_%s.pdf' % taxon
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





for taxon in ['B','C', 'D','F','J','P']:
    plot_mutation_trajectory_taxon(taxon)
