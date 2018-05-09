from __future__ import division
import os, math, numbers, itertools
import pandas as pd
import numpy as np
from string import maketrans
from collections import Counter
import  matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import sklearn.metrics.pairwise as pairwise
#import skbio.stats.ordination as ordination
#import skbio.stats.distance as distance
from matplotlib.patches import Polygon
import matplotlib as mpl
import pylab as P
import seaborn as sns

from statsmodels.graphics.factorplots import interaction_plot


mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")

strain_colors = {'Janthinobacterium':'indigo', 'Caulobacter':'lightblue', \
         'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'cyan',
         'Pedobacter': 'darkred', 'Bacillus_spoA':'darkblue'}

species_dict = {'B': 'Bacillus', 'C':'Caulobacter', 'D':'Deinococcus', \
    'F': 'Pedobacter', 'J': 'Janthinobacterium', 'P':'Pseudomonas'}

treatment_dict = {'0': '1-day', '1': '10-day', '2': '100-day'}

def CV_KDE(oneD_array):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(kernel='exponential'),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple


def AFS():
    #strains = ['B', 'C', 'D', 'F', 'J', 'P', 'S']
    strains = ['S']
    colors = {'1':'cyan', '2':'lightblue', \
             '3': 'red', '4':'darkgreen', '5':'indigo'}
    for strain in strains:
        path = mydir + 'data/breseq_output_gbk_essentials_split_clean_merged_unique_merged/Strain_' + strain + '_SNP.txt'
        if os.path.exists(path) != True:
            continue
        IN = pd.read_csv(path, sep = '\t', header = 'infer', low_memory=False)
        out_path = mydir + 'figs/SFS/' + strain
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        samples = [x for x in IN.columns if 'frequency_L' in x]
        days = set([x.split('_')[-1] for x in samples])
        for day in list(days):
            out_path_day = mydir + 'figs/SFS/' + strain + '/' + day
            if not os.path.exists(out_path_day):
                os.makedirs(out_path_day)
        for sample in samples:
            sfs = IN[sample].values
            sfs = sfs[~np.isnan(sfs)]
            sfs = sfs[(sfs != float(1))]
            if len(sfs) == 0:
                continue
            print sample
            #fig = plt.figure()
            #plt.hist(sfs, bins=30, alpha = 0.8,  normed = True)
            #plt.title(sample + ' dist. of coverage')
            #plt.xlabel('Site frequency', fontsize=14)
            #plt.ylabel('Probability', fontsize=14)
            #plt.xlim([0, max(120,  x_mean  + (x_mean * 2)) ])
            #fig.tight_layout()
            out_plot = mydir + 'figs/SFS/' + strain + '/' + sample.split('_')[-1] + '/' + sample + '.png'
            #fig.savefig(out_plot, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
            #plt.close()

            fig, ax = plt.subplots()
            #weights = np.ones_like(sfs)/len(sfs)
            #ax.hist(sfs, 50, fc=colors[rep_number], histtype='stepfilled',
            #    label='Replicate ' + rep_number, alpha=0.5, weights= weights)
            #ax.hist(sfs, 50, histtype='stepfilled', alpha=0.5, weights= weights)
            ax.hist(sfs, 40, histtype='stepfilled', alpha=0.8, normed= True)
            #print sfs
            plt.xlim([0.01,0.99])
            #ax.legend()
            ax.set_xlabel('Site frequency', fontsize = 16)
            ax.set_ylabel('Fraction of sites', fontsize = 16)
            #title = 'Site frequency spectra for ' + treatment_dict[rep[-3]] \
            #+ ' ' + species_dict[rep[-2]]
            #fig.suptitle(title, fontsize=15)
            fig.tight_layout()
            fig.savefig(out_plot, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
            plt.close()

    #KDE = CV_KDE(afs)
    #fig, ax = plt.subplots()
    #ax.plot(KDE[0], KDE[1], linewidth=3, alpha=0.5, label='bw=%.2f' % KDE[2])
    #weights1 = np.ones_like(afs1)/len(afs1)
    #weights2 = np.ones_like(afs2)/len(afs2)
    #ax.hist(afs1, 50, fc='b', histtype='stepfilled', alpha=0.5)
    #ax.hist(afs2, 50, fc='r', histtype='stepfilled', alpha=0.5)
    #plt.xlim([0.01,0.99])
    #ax.set_xlabel('Site frequency', fontsize = 16  )
    #ax.set_ylabel('Number of sites', fontsize = 16)
    #title = 'Site frequency spectreum of '
    #ax.text(2, 1650, 'Site-specific substitutions', fontsize=15)
    #plt.savefig(mydir + 'figs/SFS/D100/D/test.png', dpi=600)
    #plt.close()

    #IN[np.isfinite(IN['frequency_L0D5'])]


#def SFS_figs():

def K_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    K_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        K_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    #K_values = np.log10(K_values)
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], K_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_yscale("log", nonposy='clip')
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Substitutions', fontsize = 20)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    #ax.set_xticklabels(treatments, fontsize = 8)
    #ax.get_xaxis().tick_bottom()
    #ax.get_yaxis().tick_left()
    #ax.set_ylim(1.15, 1.7)
    #ax.set_title('Genera')
    #ax.set_ylabel('Substitutions')
    #plt.plot([1, 1, 7, 7] , [1.6, 1.65, 1.65, 1.6], lw=1.5, c='k')
    #plt.text(4, 1.66, "****", ha='center', va='bottom', color='k', fontsize = 17)
    fig.savefig(mydir + 'figs/K.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def poly_fig(variable):
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    poly_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            poly_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)][variable].tolist() )
            #treat_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)]['treatment'].tolist() )
            #strain_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)]['strain'].tolist() )
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])
    poly_values_per_gen = []
    for i, poly_value in enumerate(poly_values):
        if variable == 'W_T_L':
            poly_value = [x for x in poly_value if x != 0]
        poly_values_per_gen.append(poly_value)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(poly_values_per_gen, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    if variable == 'W_T_L':
        ax1.set_ylabel('Segregating sites, ' +  r'$\theta_{W}$' +  \
        ' \n per base', fontsize = 20)
    elif variable == 'pi_L':
        ax1.set_ylabel('Nucleotide diversity, ' +  r'$\pi$' +  \
        '\n per base', fontsize = 20)

    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*6
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(poly_values_per_gen[i])],
                 color='w', marker='*', markeredgecolor='k')
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.820, 'Pseudomonas',
                backgroundcolor='darkgreen',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.865, 'Pedobacter',
                backgroundcolor='darkred',
                color='white', weight='roman', size='x-small')

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)


    ax1.set_yscale("log", nonposy='clip')
    fig.savefig(mydir + 'figs/' + variable +'.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def pi_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], pi_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Nucleotide diversity (pi)', fontsize = 20)
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/pi.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def W_theta_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    W_theta_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        W_theta_values.append( IN[IN.strain == strain]['W_T_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], W_theta_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Wattersons theta', fontsize = 20)
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/WT.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def T_D_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', \
        sep = ' ', header = 'infer')
    #IN = IN.sort(['treatment', 'strain'], ascending=[True, True])
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    T_D_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            T_D_value =  IN[(IN.strain == strain) & (IN.treatment == treatment)]['T_D'].tolist()
            T_D_values.append([x for x in T_D_value if np.isnan(x) == False])
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])


    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = plt.boxplot(T_D_values, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    ax1.set_ylabel(r'$D_{T }$', fontsize = 26)

    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*6
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(T_D_values[i])],
                 color='w', marker='*', markeredgecolor='k')
    #labels = [item.get_text() for item in ax1.get_xticklabels()]
    #labels_new = []
    #for i, label in enumerate(labels):
    #    if i % 3 == 0:
    #        labels_new.append('1-Day')
    #    elif i % 3 == 1:
    #        labels_new.append('10-Day')
    #    elif i % 3 == 2:
    #        labels_new.append('100-Day')

    #ax1.set_xticklabels(labels_new, fontsize = 18)

    #for label in ax1.get_xmajorticklabels():
    #    label.set_rotation(60)
    #    label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
    #    label.set_horizontalalignment("right")

    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.820, 'Pseudomonas',
                backgroundcolor='darkgreen',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.865, 'Pedobacter',
                backgroundcolor='darkred',
                color='white', weight='roman', size='x-small')


    fig.savefig(mydir + 'figs/T_D.png', bbox_inches='tight',  dpi = 600)
    plt.close()



def T_D_figsss():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', \
        sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    W_theta_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        print IN[IN.strain == strain]['T_D'].tolist()
        W_theta_values.append( IN[IN.strain == strain]['T_D'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], W_theta_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Tajimas D', fontsize = 20)
    #plt.ticklabel_format(style='sci', axis='y')
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/TD.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def mut_figsss():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/mut_bias.txt', \
        sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    #treatments = list(IN.treatment.unique())
    #W_theta_values = []
    #treat_values = []
    #strain_values = []
    #colors = []
    for strain in strains:
        treatments = list(IN.treatment.unique())
        W_theta_values = []
        treat_values = []
        strain_values = []
        colors = []
        W_theta_values.append( IN[IN.strain == strain]['m_sample_ma'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
        fig, ax = plt.subplots()
        ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
        for i, color in enumerate(colors):
            ax.plot(treat_values[i], W_theta_values[i], marker='o', linestyle='', \
                ms=14, label=strain, color = color, alpha = 0.9)

        #ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
        fig.canvas.draw()
        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels[1] = '1-Day'
        labels[5] = '10-Day'
        labels[9] = '100-Day'
        plt.axhline(y=1, color='grey', linestyle='--', lw = 4)
        plt.title(strain, fontsize = 22)
        plt.ylim([0,10])
        ax.set_xticklabels(labels, fontsize = 18)
        #ax.set_ylabel(r'$m_{pop} / m_{mut},\; log_{10}$', fontsize = 20)
        ax.set_ylabel(r'$m_{pop} / m_{mut}$', fontsize = 20)
        #plt.ticklabel_format(style='sci', axis='y')
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
        fig.savefig(mydir + 'figs/mut_' + strain + '.png', bbox_inches='tight',  dpi = 600)
        plt.close()


def mut():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/mut_bias.txt', \
        sep = ' ', header = 'infer')
    #IN = IN.sort(['treatment', 'strain'], ascending=[True, True])
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    T_D_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            T_D_value =  IN[(IN.strain == strain) & (IN.treatment == treatment)]['m_sample_ma'].tolist()
            print strain
            print T_D_value
            T_D_values.append([x for x in T_D_value if np.isnan(x) == False])
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = plt.boxplot(T_D_values, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    ax1.set_ylabel(r'$D_{T }$', fontsize = 26)
    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*4
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(3):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(T_D_values[i])],
                 color='w', marker='*', markeredgecolor='k')
    #labels = [item.get_text() for item in ax1.get_xticklabels()]
    #labels_new = []
    #for i, label in enumerate(labels):
    #    if i % 3 == 0:
    #        labels_new.append('1-Day')
    #    elif i % 3 == 1:
    #        labels_new.append('10-Day')
    #    elif i % 3 == 2:
    #        labels_new.append('100-Day')

    #ax1.set_xticklabels(labels_new, fontsize = 18)

    #for label in ax1.get_xmajorticklabels():
    #    label.set_rotation(60)
    #    label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
    #    label.set_horizontalalignment("right")

    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)
    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    fig.savefig(mydir + 'figs/mut.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def mut_B_S():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/mut_bias.txt', \
        sep = ' ', header = 'infer')
    #IN = IN.sort(['treatment', 'strain'], ascending=[True, True])
    strains = ['Bacillus', 'Bacillus_spoA']
    treatments = list(IN.treatment.unique())
    ma_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            ma_value =  IN[(IN.strain == strain) & (IN.treatment == treatment)]['m_sample_ma'].tolist()
            ma_values.append([x for x in ma_value if np.isnan(x) == False])
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    print strain_values
    print ma_values
    bp = plt.boxplot(ma_values, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    ax1.set_ylabel(r'$D_{T}$', fontsize = 26)
    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*2
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(3):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(ma_values[i])],
                 color='w', marker='*', markeredgecolor='k')
    #labels = [item.get_text() for item in ax1.get_xticklabels()]
    #labels_new = []
    #for i, label in enumerate(labels):
    #    if i % 3 == 0:
    #        labels_new.append('1-Day')
    #    elif i % 3 == 1:
    #        labels_new.append('10-Day')
    #    elif i % 3 == 2:
    #        labels_new.append('100-Day')

    #ax1.set_xticklabels(labels_new, fontsize = 18)

    #for label in ax1.get_xmajorticklabels():
    #    label.set_rotation(60)
    #    label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
    #    label.set_horizontalalignment("right")

    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)
    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    fig.savefig(mydir + 'figs/mut_B_S.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def pi_vs_k2():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    k_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
        k_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        #ax.scatter(k_values[i], pi_values[i], marker='o', linestyle='', \
        #    ms=6, label=strains[i], color = color, alpha = 0.9)
        print k_values[i]
        print pi_values[i]
        ax.scatter(k_values[i], pi_values[i], marker='o', \
            label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    #labels = [item.get_text() for item in ax.get_xticklabels()]
    #labels[1] = '1-Day'
    #labels[3] = '10-Day'
    #labels[5] = '100-Day'
    ax.set_ylabel('Substitutions per-site (K)', fontsize = 20)
    ax.set_ylabel('Nucleotide diversity per-site (pi)', fontsize = 20)
    #ax.set_xscale("log", nonposy='clip')
    #ax.set_yscale("log", nonposy='clip')
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/pi_v_k.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def pi_vs_k():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    k_values = []
    treat_values = []
    strain_values = []
    colors = []
    print IN[IN.strain == 'Pedobacter']
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
        k_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    for i, color in enumerate(colors):
        ax.plot(k_values[i], pi_values[i], marker='o', alpha = 0.8, \
            linestyle='', ms=12, c = color)
        #ax.scatter(k_values[i], pi_values[i], marker='o', \
        #    label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    #print "slope = " + str(slope)
    #print "r2 = " + str(r_value**2)
    #print "p = " + str(p_value)
    #predict_y = intercept + slope * x
    #pred_error = y - predict_y
    #degrees_of_freedom = len(x) - 2
    #residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    ax.set_xlabel('Substitutions per-site (K)', fontsize = 20)
    ax.set_ylabel('Nucleotide diversity per-site (pi)', fontsize = 20)
    ax.set_xscale("log", nonposy='clip')
    sax.set_yscale("log", nonposy='clip')

    fig.savefig(mydir + 'figs/pi_v_k.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def sample_by_gene_dissimilarity(day, strain):
    path = mydir + 'data/gene_by_sample/' + strain + '/' + \
        day + '/sample_by_gene_multiple_hits.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer', index_col = 0)
    if day == 'D100' and strain == 'B':
        IN = IN.drop(['frequency_L2B3', 'frequency_L1B1'], axis=0)
    IN = IN[(IN.T != 0).any()]
    rows = IN.index
    IN_np = IN.as_matrix()
    distance_sklearn = pairwise.pairwise_distances(IN_np, metric='braycurtis')
    distance_skbio = distance.DistanceMatrix(distance_sklearn, ids=rows)
    pcoa = ordination.PCoA(distance_skbio)
    pcoa_results = pcoa.scores()
    #print pcoa_results.__dict__
    names = [i.split('_')[1] for i in pcoa_results.site_ids]
    x = pcoa_results.site[:,0]
    y = pcoa_results.site[:,1]
    zipped = zip(x, y, names)
    fig, ax = plt.subplots()
    treatments = [0, 1, 2]
    treatment_names = ['1-Day', '10-Day', '100-Day']
    treatment_colors = ['#87CEEB', '#FFA500', '#FF6347']
    for treatment in treatments:
        zipped_treatment = [k for k in zipped if str(treatment) in k[2]]
        x_treatment = [k[0] for k in zipped_treatment]
        y_treatment = [k[1] for k in zipped_treatment]
        ax.plot(x_treatment, y_treatment, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=treatment_names[treatment], \
            c = treatment_colors[treatment])
    for item in zipped:
        ax.annotate(item[2], (item[0],item[1]))
        round(14.22222223, 2)
    xlabel = 'PCoA 1 (' +  str(round(pcoa_results.proportion_explained[0] * 100, 1))  + '%)'
    ylabel = 'PCoA 2 (' +  str(round(pcoa_results.proportion_explained[1] * 100, 1))  + '%)'
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    title = 'Day 100 ' + species_dict[strain]
    plt.title(title, fontsize = 22)
    fig.savefig(mydir + 'figs/pcoa/D100_' + strain + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def GMD_hist(day, strain):
    path = mydir + 'data/gene_by_sample/' + strain + '/' + day + \
        '/sample_by_gene_multiple_hits.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer', index_col = 0)
    sample_names = IN.index.values
    day_1 = [x for x in sample_names if x[-3] == str(0) ]
    day_10 = [x for x in sample_names if x[-3] == str(1) ]
    day_100 = [x for x in sample_names if x[-3] == str(2) ]
    day_1_sum = IN.loc[day_1,:].sum(axis=0)
    day_10_sum = IN.loc[day_10,:].sum(axis=0)
    day_100_sum = IN.loc[day_100,:].sum(axis=0)
    day_1_sum.name = '1-Day'
    day_10_sum.name = '10-Day'
    day_100_sum.name = '100-Day'

    merged = pd.concat([day_1_sum, day_10_sum, day_100_sum], axis=1)
    idx = merged.sum(axis=1).sort_values(ascending=False).head(30).index
    merged_sort_descending = merged.ix[idx]
    gene_names = merged_sort_descending.index.values

    fig, ax = plt.subplots()
    bar_width = 0.65
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(merged_sort_descending['100-Day']))]
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 100-day data
            merged_sort_descending['100-Day'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='100-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#FF6347')

    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 10-day data
            merged_sort_descending['10-Day'],
            # set the width
            width=bar_width,
            # with 100-day on the bottom
            bottom=merged_sort_descending['100-Day'],
            # with the label mid score
            label='10-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#FFA500')

    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 1-day data
            merged_sort_descending['1-Day'],
            # set the width
            width=bar_width,
            # with pre_score and mid_score on the bottom
            bottom=[i+j for i,j in zip(merged_sort_descending['100-Day'],\
                merged_sort_descending['10-Day'])],
            # with the label post score
            label='1-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#87CEEB')

    # set the x ticks with names
    plt.xticks(tick_pos, gene_names)
    ax.set_ylabel("Number of mutations", fontsize = 16 )
    #ax.set_xlabel("Gene names")
    plt.legend(loc='upper right')
    for label in ax.get_xmajorticklabels():
        label.set_rotation(60)
        label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
        label.set_horizontalalignment("right")
    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)

    # Set a buffer around the edge
    plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

    title = 'Day ' + day[1:] + ' ' + species_dict[strain]
    plt.title(title, fontsize = 22)
    fig_dir = mydir + 'figs/GMD/' + day
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_out = fig_dir + '/' + strain + '.png'
    #plt.savefig(fig_out, dpi=600)
    fig.savefig(fig_out, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def get_coverage_clean_figs(sigmas = 3):
    days = ['D100', 'D200', 'D300']
    out_summary = open(mydir + 'data/coverage_clean/coverage_clean_summary.txt', 'w')
    headers = ['line', 'cov_mean', 'cov_std', 'percent_below_50', 'resequence']
    out_summary.write('\t'.join(headers) + '\n')
    #print>> out_summary, 'line', 'cov_mean', 'cov_std', 'percent_below_100', 'resequence'
    for day in days:
        in_dir = mydir + 'data/coverage_clean/' + day + '/'
        out_dir = mydir + 'figs/coverage_clean_figs/' + day + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for filename in os.listdir(in_dir):
            if filename.endswith(".txt"):
                sample = filename.split('_')[1]
                print sample
                cov_df = pd.read_csv(os.path.join(in_dir, filename), sep = ' ')
                x = cov_df.cov_mean.values
                x_mean = np.mean(x)
                if len(x) == 0:
                    continue
                x_under_cov = [x_i for x_i in x if x_i < 50]
                alpha = (len(x_under_cov) + 1) / len(x)

                fig = plt.figure()
                plt.hist(x, bins=200, alpha = 0.8, normed=True)
                plt.title(sample + ' dist. of coverage')
                plt.xlabel('Coverage', fontsize=14)
                plt.ylabel('Probability', fontsize=14)
                plt.xlim([0, max(120,  x_mean  + (x_mean * 2)) ])
                plt.axvline(x=100, c = 'grey', linestyle = '--', lw = 3)
                plt.axvline(x=x_mean, c = 'black', linestyle = '-', lw = 3)
                # std deviation in poisson process = sqrt(lambda)
                plt.axvline(x=x_mean + (np.sqrt(x_mean) * 2), c = 'black', linestyle = ':', lw = 3)
                plt.axvline(x=x_mean - (np.sqrt(x_mean) * 2), c = 'black', linestyle = ':', lw = 3)
                fig.tight_layout()
                fig.savefig(out_dir + sample + '.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
                plt.close()

                #if x_mean -  (np.sqrt(x_mean) * 2) > 100:
                #    reseq = False
                #else:
                #    reseq = True
                #print alpha
                if alpha < 0.25:
                    reseq = False
                else:
                    reseq = True
                #print>> out_summary, sample, str(x_mean), str(np.std(x)), str(alpha), str(reseq)
                data = [sample, str(x_mean), str(np.std(x)), str(alpha), str(reseq)]
                out_summary.write('\t'.join(data) + '\n')

    out_summary.close()

class make_muller_plots:

    def __init__(self, strain):
        self.strain = strain

    def stacked_trajectory_plot(self, name = 'stacked_trajectory', xlabel="generation"):
        IN = pd.read_csv(mydir + 'data/breseq_output_gbk_essentials_split_clean_merged_unique_merged/' + \
            'Strain_' + self.strain + '_SNP.txt',  sep = '\t', low_memory=False)
        colors_lighter = ["#A567AF", "#8F69C1", "#8474D1", "#7F85DB", "#7F97DF",\
                        "#82A8DD", "#88B5D5", "#8FC0C9", "#97C8BC", "#A1CDAD", \
                        "#ACD1A0", "#B9D395", "#C6D38C", "#D3D285", "#DECE81", \
                        "#E8C77D", "#EDBB7A", "#EEAB77", "#ED9773", "#EA816F", \
                        "#E76B6B"]
        pops = [x.split('_')[1] for x in IN.columns if 'frequency_' in x]
        for pop in list(set(pops)):
            if pop != 'frequency':
                time_points = [x for x in IN.columns if 'frequency_' + pop in x]
                time_points_x = [int(x.split('_')[2][1:]) for x in time_points]
                time_points_x.insert(0, 0)
                IN_time_points = IN[time_points]
                #trajectories = IN_time_points.dropna(thresh=2)
                trajectories = IN_time_points.dropna(how='all')
                #print trajectories.loc[trajectories['frequency_L0P2_D200'] > 0.9]
                # NaNs should be counted as 0, since we have fixed mutations in the df
                #trajectories['D0'] = 0
                trajectories.insert(0, 'D0', 0)
                #print trajectories
                trajectories = trajectories.fillna(0)
                trajectories = trajectories.values
                if trajectories.shape[0] == 0:
                    continue

                fig = plt.figure()
                for i in range(trajectories.shape[0]):
                    plt.plot(time_points_x, trajectories[i,], linestyle='-', \
                        marker='o', alpha = 0.6)
                print pop
                plt.ylim(0, 1)
                plt.xlim(0, 300)
                plt.ylabel("Frequency")
                plt.xlabel('Days')
                plt.title('Population ' + pop)
                plt.xticks(np.arange(min(time_points_x), max(time_points_x)+1, 100))
                out_dir = mydir + 'figs/muller_plots/' + self.strain
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                plt.savefig(out_dir + '/' + pop + '.png', bbox_inches='tight',  dpi = 600)
                plt.close()


def euc_dist(strain):
    IN = pd.read_csv(mydir + 'data/euclidean_distance/' +  strain + '.txt', sep = '\t', low_memory=False)
    IN_0_100 =  IN[IN['Line'].str.contains("L0B") & (IN['Time1'] == 100)]
    IN_1_100 =  IN[IN['Line'].str.contains("L1B") & (IN['Time1'] == 100)]
    IN_2_100 =  IN[IN['Line'].str.contains("L2B") & (IN['Time1'] == 100)]
    print np.mean(IN_0_100.Distance.values), np.std(IN_0_100.Distance.values)
    print np.mean(IN_1_100.Distance.values), np.std(IN_1_100.Distance.values)
    print np.mean(IN_2_100.Distance.values), np.std(IN_2_100.Distance.values)


def pi_b_s(param = 'pi_L'):
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    IN = IN[IN.pi_L != 0]
    b_1 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 0)]
    b_10 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 1)]
    b_100 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 2)]
    s_1 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 0)]
    s_10 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 1)]
    s_100 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 2)]
    #mean_100_1day = mean_100.loc[mean_100['TransferTime'] == 1]
    #mean_100_10day = mean_100.loc[mean_100['TransferTime'] == 10]

    data_to_plot = [b_1[param], b_10[param], b_100[param], s_1[param], s_10[param], s_100[param]]
    data_to_plot = [np.log10(x) for x in data_to_plot]
    print IN
    #print data_to_plot[1]#data_to_plot = [x[(x.T != 0).any()]]
    print data_to_plot
    #data_to_plot = [np.log10(x) for x in data_to_plot ]
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))
    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot)
    ## add patch_artist=True option to ax.boxplot()
    ## to get fill color
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    ## Custom x-axis labels
    #ax.set_xticklabels(['1-Day', '10-Day'])
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylabel(r'$\pi$', fontsize = 25)

    out_dir = mydir + 'figs/pi_b_s.png'
    plt.savefig(out_dir, bbox_inches='tight',  dpi = 600)
    plt.close()
    #(df["B"] > 50) & (df["C"] == 900)

def td_b_s(param = 'T_D'):
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    IN = IN[IN.pi_L != 0]
    b_1 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 0)]
    b_10 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 1)]
    b_100 = IN.loc[(IN['strain'] == 'Bacillus') & (IN['treatment'] == 2)]
    s_1 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 0)]
    s_10 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 1)]
    s_100 = IN.loc[(IN['strain'] == 'Bacillus_spoA') & (IN['treatment'] == 2)]
    #mean_100_1day = mean_100.loc[mean_100['TransferTime'] == 1]
    #mean_100_10day = mean_100.loc[mean_100['TransferTime'] == 10]

    data_to_plot = [b_1[param], s_1[param], b_10[param], s_10[param], b_100[param], s_100[param]]
    #print data_to_plot[1]#data_to_plot = [x[(x.T != 0).any()]]
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))
    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot)
    ## add patch_artist=True option to ax.boxplot()
    ## to get fill color
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    ## Custom x-axis labels
    #ax.set_xticklabels(['1-Day', '10-Day'])
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    #ax.set_yscale("log", nonposy='clip')
    ax.set_ylabel(r'$\pi$', fontsize = 25)

    out_dir = mydir + 'figs/td_b_s.png'
    plt.savefig(out_dir, bbox_inches='tight',  dpi = 600)
    plt.close()


def beta_fig():
    IN = pd.read_csv(mydir + 'data/betas.txt', sep = '\t', header = 'infer')
    IN['Response'] = IN.sum(axis = 1).values
    print IN

    B0 = IN.index[IN.index.str.contains('L0B')].values
    B1 = IN.index[IN.index.str.contains('L1B')].values
    B2 = IN.index[IN.index.str.contains('L2B')].values
    S0 = IN.index[IN.index.str.contains('L0S')].values
    S1 = IN.index[IN.index.str.contains('L1S')].values
    S2 = IN.index[IN.index.str.contains('L2S')].values
    #betas = np.append(B0, B1, B2, S0, S1, S2)
    #betas = np.hstack(( B0, B1, B2, S0, S1, S2 )).ravel()

    strain = np.asarray((['B'] * (len(B0) + len(B1) + len(B2) ))  + (['S'] * (len(S0) + len(S1) + len(S2) )))
    time = np.asarray( ([1] * len(B0))  + ([10] * len(B1)) + ([100] * len(B2)) + ([1] * len(S0))  + ([10] * len(S1)) + ([100] * len(S2)) )
    IN['Strain'] = strain
    IN['Time'] = time
    #df = pd.DataFrame({'Betas':betas, 'Strain':strain, 'Time':time})

    sns_plot = sns.factorplot(kind='box',        # Boxplot
               y='Response',       # Y-axis - values for boxplot
               x='Time',        # X-axis - first factor
               hue='Strain',         # Second factor denoted by color
               data=IN,        # Dataframe
               size=8,            # Figure size (x100px)
               aspect=1.5,        # Width = size * aspect
               legend_out=False,
               ci = "sd")
    #sns_plot.set(xlabel='Time', ylabel=r'$\beta_{Time:Strain}$', font= 16)
    sns_plot.savefig(mydir + "figs/betas.png")



def beta_fig_inter():
    IN = pd.read_csv(mydir + 'data/betas.txt', sep = '\t', header = 'infer')
    IN['Response'] = IN.sum(axis = 1).values

    B0 = IN.index[IN.index.str.contains('L0B')].values
    B1 = IN.index[IN.index.str.contains('L1B')].values
    B2 = IN.index[IN.index.str.contains('L2B')].values
    S0 = IN.index[IN.index.str.contains('L0S')].values
    S1 = IN.index[IN.index.str.contains('L1S')].values
    S2 = IN.index[IN.index.str.contains('L2S')].values

    strain = np.asarray((['B'] * (len(B0) + len(B1) + len(B2) ))  + (['S'] * (len(S0) + len(S1) + len(S2) )))
    time = np.asarray( ([1] * len(B0))  + ([10] * len(B1)) + ([100] * len(B2)) + ([1] * len(S0))  + ([10] * len(S1)) + ([100] * len(S2)) )
    IN['Strain'] = strain
    IN['Time'] = time
    #print response
    response = IN['Response'].values


    fig, ax = plt.subplots(figsize=(9, 6))
    fig = interaction_plot(x=time, trace=strain, response=response,
                           colors=['#87CEEB', '#FF6347'], markers=['D', '^'], ms=10, ax=ax)
    #fig = int_plot( x=time, trace=strain, response=response, errorbars=True, errorbartyp='std',  ax=ax)
    L=plt.legend()
    L.get_texts()[0].set_text('wt')
    L.get_texts()[1].set_text(r'$\Delta spo0A$')
    plt.title('PERMANOVA interaction plot', fontsize = 20)
    plt.xlabel('Transfer time (days), ' + r'$log_{10}$', fontsize = 18)
    plt.ylabel('Mean dissimilarity', fontsize = 18)
    ax.set_xscale("log")#, nonposy='clip')
    fig_name = mydir + 'figs/RegressInterPlot.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

    #x_M = IN.loc[IN['Gender'] == 'M'].Midparent
    #x_F = IN.loc[IN['Gender'] == 'F'].Midparent
    #y_M = IN.loc[IN['Gender'] == 'M'].Height
    #y_F = IN.loc[IN['Gender'] == 'F'].Height

    #fig = plt.figure()
    #plt.scatter(x_M, y_M, c='#87CEEB', marker='o', label='Men')
    #plt.scatter(x_F, y_F, c='#FF6347', marker='o', label='Women')
    #y_pred_F = mod1.params[0] + mod1.params[1] * 0 + mod1.params[2] * midparent
    #y_pred_M = mod1.params[0] + mod1.params[1] * 1 + mod1.params[2] * midparent
    #plt.plot(midparent, y_pred_F, 'k-', lw = 5, c = 'black', label = '_nolegend_' )
    #plt.plot(midparent, y_pred_F, 'k-', lw = 2, c = '#FF6347', label = '_nolegend_')
    #plt.plot(midparent, y_pred_M, 'k-', lw = 5, c = 'black', label = '_nolegend_' )
    #plt.plot(midparent, y_pred_M, 'k-', lw = 2, c = '#87CEEB', label = '_nolegend_')

    #plt.plot(midparent, y_pred_F, c = '#FF6347')
    #plt.plot(midparent, y_pred_M, c = '#87CEEB')

    #fig_name = mydir + 'Figures/galtonRegressInter.png'
    #fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.close()



#poly_fig()
#K_fig()
#pi_fig()
#W_theta_fig()
#T_D_fig()
#pi_vs_k()
#poly_fig('pi_L')
#poly_fig('W_T_L')
#poly_fig('k_L')
#mut()
#mut_figsss()
beta_fig_inter()

#get_coverage_clean_figs()
#AFS()
#strains = ['B', 'C', 'D', 'F', 'J', 'P']
#for strain in strains:
#    make_muller_plots(strain).stacked_trajectory_plot()
#make_muller_plots('P').stacked_trajectory_plot()
#euc_dist('B')
#mut_B_S()
#beta_fig()
