from __future__ import division
import os, math
import pandas as pd
import numpy as np
from mpmath import nsum, inf

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")

def calculate_null_dist():
    return(232)


def heavside_step_fxn(x):
    if x >= 0:
        return(1)
    else:
        return(0)



def inner_numm_m_dist(n, m, n_tot, n_i, L_mean, L_i, N_genes):
    heavside_x = (n_i * (L_mean/L_i)) - m
    term1 = (n / n_tot) * heavside_step_fxn(heavside_x)
    term2 = (((n_tot*L_i) / (L_mean*N_genes)) ** n) / math.factorial(n)
    term3 = math.exp(- ((n_tot*L_i) / (L_mean*N_genes)) )
    print(term1 * term2 * term3)
    return(term1 * term2 * term3)


def get_null():
    #m_path = mydir + 'data/gene_by_sample/B_S/gene_by_sample_m.txt'
    count_path = mydir + 'data/gene_by_sample/B_S/sample_by_gene.txt'
    df = pd.read_csv(count_path, sep = '\t', header = 'infer', index_col = 0)
    df = df[df.index.to_series().str.contains('D100')]
    gene_path = mydir + "data/reference_assemblies_task2/reference_assemblies_task2_table/B.txt"
    df_gene = pd.read_csv(gene_path, sep = ' ', header = 'infer', index_col = 0)

    terms = []
    #print(df_gene)
    # merge genes as column, using gene legnth
    
    #for x in range(0),
    #print(nsum(lambda n: inner_numm_m_dist(n, 3, 25, 10, 1000, 1000, 800), [0, inf]))


    #L0B = df[df.index.to_series().str.contains('L0B')]
    #print(df_gene)






get_null()
