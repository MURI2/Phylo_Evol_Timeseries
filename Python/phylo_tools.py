from __future__ import division
import os
from collections import Counter
import matplotlib.colors as clr
import numpy as np
import colorsys


def get_path():
    return os.path.expanduser("~/GitHub/Phylo_Evol_Timeseries")



def samples_to_remove(population):
    population_dict = {'0D1':[600],
                        '0D2':[900],
                        '0D3':[900],
                        '0D4':[900]}


    if population not in population_dict:
        return None
    else:
        return population_dict[population]



def pops_to_ignore():
    return ['1C1']

def get_treatment_name(treatment):
    treatment_dict = {'0': '1-day',
                        '1': '10-day',
                        '2': '100-day'}
    return treatment_dict[treatment]



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



def get_genome_size(taxon):
    genome_size_dict = {"B": 4299822,
                        "S": 4299822,
                        "C": 4042929,
                        "D": 3284156,
                        "F": 6337316,
                        "J": 6118925,
                        "P": 6639537}

    return genome_size_dict[taxon]



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
                "F": "data/reference_assemblies_task2/Pedobacter_sp_KBS0701/GCA_005938645.1_ASM593864v1_genomic.gbff",
                "J": "data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/GCA_005937955.1_ASM593795v1_genomic.gbff",
                "P": "data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/GCA_005938045.1_ASM593804v1_genomic.gbff"}

    return ref_dict[taxon]



def get_ref_fna_dict():

    ref_dict = {"B": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna",
                "S": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna",
                "C": "data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.fna",
                "D": "data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.fna",
                "F": "data/reference_assemblies_task2/Pedobacter_sp_KBS0701/GCA_005938645.1_ASM593864v1_genomic.fna",
                "J": "data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/GCA_005937955.1_ASM593795v1_genomic.fna",
                "P": "data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/GCA_005938045.1_ASM593804v1_genomic.fna"}
    return ref_dict



def get_ref_bresq_fna_dict():

    ref_dict = {"B": "data/reference_assemblies_task2/breseq_refs/0B1_100.fasta",
                "S": "data/reference_assemblies_task2/breseq_refs/0B1_100.fasta",
                "C": "data/reference_assemblies_task2/breseq_refs/GCA_000022005.1_ASM2200v1_genomic.fna",
                "D": "data/reference_assemblies_task2/breseq_refs/0D1_100.fasta",
                "F": "data/reference_assemblies_task2/breseq_refs/GCA_005938645.1_ASM593864v1_genomic.fna",
                "J": "data/reference_assemblies_task2/breseq_refs/GCA_005937955.1_ASM593795v1_genomic.fna",
                "P": "data/reference_assemblies_task2/breseq_refs/GCA_005938045.1_ASM593804v1_genomic.fna"}
    return ref_dict




def get_to_keep():
    return ['SNP', 'INS', 'DEL']


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

def mutations_to_exclude():
    directory = os.fsencode(get_path() + '/data/pool_pop_seq/rebreseq_annotated')
    list_muts = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.gd'):
            in_df = open(os.path.join(str(directory, 'utf-8'), filename), 'r')
            for line in in_df:
                line_split = line.strip().split()
                if line_split[0] not in get_to_keep():
                    continue
                freq = [s for s in line_split if 'frequency=' in s][0].split('=')[1]
                if float(freq) == float(1):
                    list_muts.append(line_split[3] + '_' + line_split[4])
    dict_muts = Counter(list_muts)
    return [i for i in dict_muts if dict_muts[i] >= 10]



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
