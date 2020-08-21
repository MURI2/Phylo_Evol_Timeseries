from __future__ import division
import os, sys
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    Returns
    -------
    matplotlib.patches.Ellipse
    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)






treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']

spectra_dict = {}

for taxon in taxa:

    ma_mutation_spectrum = pt.get_ma_mutation_spectrum(taxon)

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



df = pd.DataFrame.from_dict(spectra_dict)
X = df.values.transpose()

X= StandardScaler().fit_transform(X) # normalizing the features, subtract mean divide by SD

normalised_df = pd.DataFrame(X,columns=df.index, index=df.columns)
pca_ = PCA(n_components=2)
principalComponents_ = pca_.fit_transform(X)

principalComponents_df = pd.DataFrame(principalComponents_, index=df.columns, columns = ['PC1', 'PC2'])

print(principalComponents_df)


fig, ax = plt.subplots(figsize=(6, 6))

for treatment in treatments:
    PCs_treatment = principalComponents_df[principalComponents_df.index.str.contains(treatment)]

    ax.scatter(PCs_treatment.PC1.values, PCs_treatment.PC2.values, \
            c=pt.get_colors(treatment), marker = 'o', s = 70, \
            edgecolors='#244162', linewidth = 0.6, alpha = 0.5, zorder=2)#, edgecolors='none'

    confidence_ellipse(PCs_treatment.PC1.values, PCs_treatment.PC2.values, ax,
        n_std=2, edgecolor=pt.get_colors(treatment), linestyle='--', lw=3)



ax.set_xlabel('PC 1 (' + str(round(pca_.explained_variance_ratio_[0]*100,2)) + '%)' , fontsize = 13)
ax.set_ylabel('PC 2 (' + str(round(pca_.explained_variance_ratio_[1]*100,2)) + '%)' , fontsize = 13)

fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/mutation_spectra_pca.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
