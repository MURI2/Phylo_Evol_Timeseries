from __future__ import division
import os
from collections import Counter
import numpy as np
import colorsys

import matplotlib.colors as clr
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from scipy.linalg import block_diag
from sklearn.metrics.pairwise import euclidean_distances

import parse_file
#import timecourse_util

from asa159 import rcont2



def get_path():
    return os.path.expanduser("~/GitHub/Phylo_Evol_Timeseries")

taxa = ['B','C','D','F','J','P']
treatments = ['0', '1', '2']
replicates = ['1','2','3','4','5']

sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

samples_to_remove = {'1B4':[900],
                        '1B5':[1000],
                        '0C1':[900],
                        '0C5':[300],
                        '1C4':[300],
                        '0D1':[600],
                        '0D2':[900],
                        '0D3':[900],
                        '0D4':[900],
                        '1J3':[400],
                        '1J4':[500],
                        '2S2':[700],
                        '2S3':[700],
                        '2S4':[700],
                        '2S5':[700]}


#populations_to_ignore = ['1D4', '2P4', '2P5', '2F1', '2F2', '2F3', '0J2', '1J3'] # ['1C1']
populations_to_ignore = ['0F3', '2F1', '2F2', '2F3', '0J2', '1J1', '1J2', '1J3', '1J4', '1J5', '1P5'] # ['1C1']

treatment_taxa_to_ignore = ['1J']





latex_formal_dict = {  'B': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610}$',
                'S': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610} \, \Delta \mathrm{spo0A} $',
                'C': r'$\mathit{Caulobacter \, crescentus} \; \mathrm{NA1000}$',
                'D': r'$\mathit{Deinococcus \, radiodurans} \; \mathrm{BAA-816}$',
                'P': r'$\mathit{Pseudomonas \,} \; \mathrm{sp. \, KBS0710}$',
                'F': r'$\mathit{Pedobacter \,} \; \mathrm{sp. \, KBS0701}$',
                'J': r'$\mathit{Janthinobacterium \,} \; \mathrm{sp. \, KBS0711}$'
                }



latex_dict = {  'B': r'$\mathit{Bacillus\, subtilis}  $',
                'S': r'$\mathit{Bacillus\, subtilis} \, \Delta \mathrm{spo0A} $',
                'C': r'$\mathit{Caulobacter \, crescentus}$',
                'D': r'$\mathit{Deinococcus \, radiodurans}$',
                'P': r'$\mathit{Pseudomonas \,} \; \mathrm{sp. \, KBS0710}$',
                'F': r'$\mathit{Pedobacter \,} \; \mathrm{sp. \, KBS0701}$',
                'J': r'$\mathit{Janthinobacterium \,} \; \mathrm{sp. \, KBS0711}$'
                }

latex_bold_dict = {'B': r'$\mathbf{\mathit{Bacillus\, subtilis} \, \mathrm{wt} }$',
                'S': r'$\mathbf{\mathit{Bacillus\, subtilis} \, \Delta \mathrm{spo0A}} $',
                'C': r'$\mathbf{\mathit{Caulobacter \, crescentus}}$',
                'D': r'$\mathbf{\mathit{Deinococcus \, radiodurans}}$',
                'P': r'$\mathbf{\mathit{Pseudomonas \,} \; \mathrm{sp. \, KBS0710}}$',
                'F': r'$\mathbf{\mathit{Pedobacter \,} \; \mathrm{sp. \, KBS0701}}$',
                'J': r'$\mathbf{\mathit{Janthinobacterium \,} \; \mathrm{sp. \, KBS0711}}$'}



#latex_genus_dict = {  'B': r'$\mathit{Bacillus} \, \mathrm{wt} $',
#                'S': r'$\mathit{Bacillus} \, \Delta \mathrm{spo0A} $',
#                'C': r'$\mathit{Caulobacter}$',
#                'D': r'$\mathit{Deinococcus}$',
#                'P': r'$\mathit{Pseudomonas}$',
#                'F': r'$\mathit{Pedobacter}$',
#                'J': r'$\mathit{Janthinobacterium}$'
#                }


latex_genus_dict = {  'B': r'$\mathit{Bacillus} $',
                'C': r'$\mathit{Caulobacter}$',
                'D': r'$\mathit{Deinococcus}$',
                'P': r'$\mathit{Pseudomonas}$',
                'F': r'$\mathit{Pedobacter}$',
                'J': r'$\mathit{Janthinobacterium}$'
                }


#latex_genus_bold_dict = {'B': r'$\mathbf{\mathit{Bacillus} }$',
#                'S': r'$\mathbf{\mathit{Bacillus} \, \Delta \mathrm{spo0A}} $',
#                'C': r'$\mathbf{\mathit{Caulobacter}}$',
#                'D': r'$\mathbf{\mathit{Deinococcus}}$',
#                'P': r'$\mathbf{\mathit{Pseudomonas}}$',
#                'F': r'$\mathbf{\mathit{Pedobacter}}$',
#                'J': r'$\mathbf{\mathit{Janthinobacterium} }$'}


latex_genus_bold_dict = {'B': r'$\mathbf{\mathit{Bacillus} }$',
                'C': r'$\mathbf{\mathit{Caulobacter}}$',
                'D': r'$\mathbf{\mathit{Deinococcus}}$',
                'P': r'$\mathbf{\mathit{Pseudomonas}}$',
                'F': r'$\mathbf{\mathit{Pedobacter}}$',
                'J': r'$\mathbf{\mathit{Janthinobacterium} }$'}




genus_dict = {'B':'Bacillus',
            'C':'Caulobacter',
            'D':'Deinococcus',
            'F':'Pedobacter',
            'J':'Janthinobacterium',
            'P':'Pseudomonas'}


def get_p_value_latex(p_value, alpha=0.05):

    if p_value < alpha:
        return r'$\mathrm{p} < 0.05$'

    else:
        return r'$\mathrm{p} \nless 0.05$'





def get_B_S_generations(strain, treatment, day_cutoff=500):

    B_S_generation_dict = {'B': { '0':3321, '1':1171, '2': 107},
                            'S': {'0':3321, '1':544, '2':163} }

    return B_S_generation_dict[strain][treatment] * (day_cutoff/1000)


def get_ma_mutation_spectrum(taxon):

    # 10^-10/site/generation

    ma_dict = {'B': {'GC_AT': 2.80,
                    'GC_TA': 0.51,
                    'GC_CG': 0.17,
                    'AT_GC': 2.16,
                    'AT_CG': 0.44,
                    'AT_TA':0.46},

                'C': {'GC_AT': 2.11,
                    'GC_TA': 0.14,
                    'GC_CG': 0.38,
                    'AT_GC': 3.87,
                    'AT_CG': 0.61,
                    'AT_TA': 0.69},

                'D': {'GC_AT': 3.01,
                    'GC_TA': 0.63,
                    'GC_CG': 0.17,
                    'AT_GC': 3.43,
                    'AT_CG': 3.19,
                    'AT_TA': 0.74},

                'J': {'GC_AT': 1.33,
                    'GC_TA': 0.19,
                    'GC_CG': 0.15,
                    'AT_GC': 0.36,
                    'AT_CG': 0.15,
                    'AT_TA': 0}}


    ma_dict_relative = {k: v / sum(ma_dict[taxon].values()) for k, v in ma_dict[taxon].items()}

    return ma_dict_relative










def run_permanova(PC_space, N_list, iter=10000):

    N_list = np.asarray(N_list)

    F_obs = get_F_2(PC_space, N_list)

    F_permute_list = []

    for i in range(iter):

        PC_space_permute = PC_space[np.random.permutation(PC_space.shape[0]),:]
        F_permute_list.append(get_F_2(PC_space_permute, N_list))

    p = len([j for j in F_permute_list if j > F_obs]) / iter

    return F_obs, p


def get_F_2(PC_space, N_list):
    '''
    Modified F-statistic from Anderson et al., 2017 doi: 10.1111/anzs.12176
    Function assumes that the rows of the count matrix are sorted by group
    i.e., group one is first N1 rows, group two is N2, etc
    '''
    #N = N1+N2
    N = sum(N_list)
    dist_matrix = euclidean_distances(PC_space, PC_space)
    A = -(1/2) * (dist_matrix ** 2)
    I = np.identity(N)
    J_N = np.full((N, N), 1)
    G = (I - ((1/N) * J_N )) @ A @ (I - ((1/N) * J_N ))
    # n matrix list
    n_list = []
    for N_i in N_list:
        n_list.append((1/N_i) * np.full((N_i, N_i), 1))
    #n1 = (1/N1) * np.full((N1, N1), 1)
    #n2 = (1/N2) * np.full((N2, N2), 1)
    #H = block_diag(n1, n2) - ((1/N) * J_N )
    H = block_diag(*n_list) - ((1/N) * J_N )
    # indicator matrices
    # get V matrices
    V_list = []
    for i in range(len(N_list)):
        if i == 0:
            U_i = np.diag( N_list[i]*[1] + sum(N_list[i+1:])*[0])
        elif i == len(N_list) - 1:
            U_i = np.diag( sum(N_list[:i])*[0] + N_list[i]*[1] )
        else:
            U_i = np.diag( sum(N_list[:i])*[0] + N_list[i]*[1] +  sum(N_list[i+1:])*[0])

        V_i = np.trace(((I - H) @ U_i @ (I - H)) @ G ) / (N_list[i]-1)
        V_list.append(V_i)

    F_2 = np.trace(H @ G) / sum( [ (1 - (N_list[i]/N) ) *  V_list[i] for i in range(len(N_list)) ]  )



    return F_2







def get_taxon_ls(taxon):

    if taxon == 'S':
        ls = ':'
    else:
        ls ='--'

    return ls



def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])



#def samples_to_remove(population):


#    if population not in population_dict:
#        return None
#    else:
#        return population_dict[population]



#'0D1':[600],
#'0D2':[900],
#'0D3':[900],
#'0D4':[900],
#'2S2':[700],
#'2S3':[700],
#'2S4':[700],
#'2S5':[700],
#'1P5':[700],
#'2P1':[500],
#'1F5':[500],
#'1F2':[700],
#'2J4':[500]



def get_treatment_name(treatment):
    treatment_dict = {'0': '1-day',
                        '1': '10-day',
                        '2': '100-day'}
    return treatment_dict[str(treatment)]



def plot_species_marker(taxon):

    plot_species_marker_dict = {"B": "o",
                                "S": "o",
                                "C": "^",
                                "D": "D",
                                "F": "P",
                                "J": "s",
                                "P": "*"}

    return plot_species_marker_dict[taxon]


def plot_species_fillstyle(taxon):

    plot_species_fillstyle_dict = {"B": "full",
                                "S": "none",
                                "C": "full",
                                "D": "full",
                                "F": "full",
                                "J": "full",
                                "P": "full"}

    return plot_species_fillstyle_dict[taxon]


def get_colors(treatment):
    get_colors_dict = {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}
    return get_colors_dict[treatment]


#def get_scatter_edge_color(strain, treatment):



def get_scatter_facecolor(taxon, treatment):

    if taxon == 'S':
        return 'white'
    else:
        return get_colors(treatment)




def get_genome_size(taxon):
    genome_size_dict = {"B": 4299822,
                        "S": 4299822,
                        "C": 4042929,
                        "D": 3284156,
                        "F": 6337316,
                        "J": 6118925,
                        "P": 6639537}

    return genome_size_dict[taxon]





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




def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = clr.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*clr.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])




def mut_freq_colormap():
    #cmap = clr.LinearSegmentedColormap.from_list('Zissou1', ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"], N=256)
    cmap = clr.LinearSegmentedColormap.from_list('Darjeeling1', ["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"], N=256)

    # sample from cmap using uniform dist b/w 0 and 1
    u = np.random.uniform()
    rgb = '#%02x%02x%02x' % tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    #tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    # RGB six digit code
    return rgb





def get_populations(taxon):

    pop_dict = {"B": ['0B1', '0B2', '0B3', '0B4', '0B5',
                        '1B1', '1B2', '1B4', '1B4', '1B5',
                        '2B1', '2B2', '2B4', '2B4', '2B5'],
                "S": ['0S1', '0S2', '0S3', '0S4', '0S5',
                        '1S1', '1S2', '1S4', '1S4', '1S5',
                        '2S1', '2S2', '2S4', '2S4', '2S5']
                }



def get_ref_gbff_dict(taxon):

    ref_dict = {"B": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.gbff",
                "S": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.gbff",
                "C": "data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCF_000022005.1_ASM2200v1_genomic.gbff",
                "D": "data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCF_000008565.1_ASM856v1_genomic.gbff",
                "F": "data/reference_assemblies_task2/Pedobacter_sp_KBS0701/GCF_005938645.2_ASM593864v2_genomic.gbff",
                "J": "data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/GCF_005937955.2_ASM593795v2_genomic.gbff",
                "P": "data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/GCF_005938045.2_ASM593804v2_genomic.gbff"}

    return ref_dict[taxon]



def get_ref_fna_dict():

    ref_dict = {"B": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna",
                "S": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna",
                "C": "data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCF_000022005.1_ASM2200v1_genomic.fna",
                "D": "data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCF_000008565.1_ASM856v1_genomic.fna",
                "F": "data/reference_assemblies_task2/Pedobacter_sp_KBS0701/GCF_005938645.2_ASM593864v2_genomic.fna",
                "J": "data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/GCF_005937955.2_ASM593795v2_genomic.fna",
                "P": "data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/GCF_005938045.2_ASM593804v2_genomic.fna"}
    return ref_dict





def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def get_bacillus_mut_bias():
    return  1.2739

def get_bacillus_mut_rate():
    return 3.35 * (10**-10)

def get_bacillus_indel_rate():
    return 1.20 * (10**-10)



def get_ts_tv_dict():

    ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V"}

    return ts_tv_dict

def get_codon_dict():
    # translation table 11
    codon_dict = {
        "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
        "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
        "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
        "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

        "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
        "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
        "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
        "CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

        "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
        "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
        "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
        "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

        "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
        "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
        "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
        "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
        }

    return codon_dict




class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list
