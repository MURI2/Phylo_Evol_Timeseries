from __future__ import division
import os, math
from collections import Counter
import numpy as np
import pandas as pd
from skbio.diversity import beta_diversity
import matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")

def cmdscale(D):
    """
    Classical multidimensional scaling (MDS)

    Parameters
    ----------
    D : (n, n) array
        Symmetric distance matrix.

    Returns
    -------
    Y : (n, p) array
        Configuration matrix. Each column represents a dimension. Only the
        p dimensions corresponding to positive eigenvalues of B are returned.
        Note that each dimension is only determined up to an overall sign,
        corresponding to a reflection.

    e : (n,) array
        Eigenvalues of B.

    """
    # Number of points
    n = len(D)

    # Centering matrix
    H = np.eye(n) - np.ones((n, n))/n

    # YY^T
    B = -H.dot(D**2).dot(H)/2

    # Diagonalize
    evals, evecs = np.linalg.eigh(B)

    # Sort by eigenvalue in descending order
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    # Compute the coordinates using positive-eigenvalued components only
    w, = np.where(evals > 0)
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)

    return Y, evals


def get_ds(df):
    num_pops = np.shape(df)[0]
    matrix = np.zeros((num_pops, num_pops))
    def calculate_bc(site_i, site_j):
        S_i = sum(site_i)
        S_j = sum(site_j)
        zip_i_j = list(zip(site_i, site_j))
        #print([x for x in zip_i_j if (x[0] > float(0)) and (x[1] > float(0))])
        C_ij = sum([min(x) for x in zip_i_j if (x[0] > float(0)) and (x[1] > float(0))])
        BC_ij = 1 - ((2*C_ij) /  (S_i + S_j)  )
        return(BC_ij)

    for i, d_i in enumerate(df[:,1:]):
        for j, d_j in enumerate(df[:i,1:]):
            bc = calculate_bc(df[i,1:], df[j,1:])
            matrix[i,j] = bc
            matrix[j,i] = bc
    return(matrix)



def merge_B_S_sample_by_gene_matrix():
    B_path = mydir + 'data/gene_by_sample/B/sample_by_gene.txt'
    S_path = mydir + 'data/gene_by_sample/S/sample_by_gene.txt'
    B = pd.read_csv(B_path, sep = '\t', header = 'infer', index_col = 0)
    S = pd.read_csv(S_path, sep = '\t', header = 'infer', index_col = 0)
    B_S_merged = B.append(S)
    B_S_merged = B_S_merged.fillna(0.0)
    # add columns of zero for genes with no observed mutations
    gene_path = mydir + 'data/reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    IN_gene = pd.read_csv(gene_path, sep = ' ', header = 'infer')
    to_add = set(IN_gene.LocusTag.values).difference(set(B_S_merged.columns.values))
    for gene in to_add:
        B_S_merged[gene] = float(0)
    # use set theory to get genes that aren't listed in matrix
    # add those genes as empty columns
    OUTname = mydir + 'data/gene_by_sample/B_S/sample_by_gene.txt'
    B_S_merged.to_csv(OUTname, sep = '\t', index = True)


def get_sample_by_gene_matrix_gscore(count_matrix):
    gene_path = mydir + 'data/reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    IN_gene = pd.read_csv(gene_path, sep = ' ', header = 'infer')
    # remove intergenic regions
    cols = [c for c in count_matrix.columns if '/' not in c]
    count_matrix = count_matrix[cols]
    L_tot = sum(IN_gene.Size.values)
    size_dict = dict(zip(IN_gene.LocusTag, IN_gene.Size))
    rel_sizes = np.asarray([size_dict[x] for x in count_matrix.columns]) / L_tot
    for index, row in count_matrix.iterrows():
        E_i = int(sum(row)) * rel_sizes
        E_i_and_G_i = list(zip(E_i, row))
        G_i = [2*i[1]*np.log(i[1]/i[0]) for i in E_i_and_G_i]
        G_i = [0 if math.isnan(y) else y for y in G_i]

        # make negatives zero for now. Need to worry about this later
        #num_neg = [x for x in G_i if x < 0]
        G_i = [0 if z < 0 else z for z in G_i]
        #value_when_true if condition else value_when_false
        count_matrix.loc[index,:] = G_i

    #OUTname = mydir + 'gene_by_sample/' + strain + '/sample_by_gene_Gscore.txt'
    #IN_gbp.to_csv(OUTname, sep = '\t', index = True)
    #print(count_matrix['frequency_L1B3_D200'])
    #print(count_matrix)
    return(count_matrix)


def get_pcoa(df):
    # remove columns containing only zeros
    df_no0 = df.loc[:, (df != 0).any(axis=0)]
    # only keep pops from day 100
    ids = df_no0.index.values
    data = df_no0.values
    ds = get_ds(data)
    pcoa = cmdscale(ds)
    Y = pd.DataFrame(pcoa[0])
    Y['pops'] = ids
    Y = Y.set_index('pops')
    return([Y, pcoa[1]])


def get_random_gene_by_pop(matrix, gene_table):
    genes = matrix.columns.values
    index_map = {v: i for i, v in enumerate(genes)}
    L_tot = sum(gene_table.Size.values)
    rel_sizes = dict(zip(gene_table.LocusTag, gene_table.Size / L_tot))
    rel_sizes_sorted = (sorted(rel_sizes.items(), key=lambda pair: index_map[pair[0]]))
    genes = [x[0] for x in rel_sizes_sorted]
    pmf = [x[1] for x in rel_sizes_sorted]
    rndm_dict = {}
    for gene in genes:
        rndm_dict[gene] = {}
        for pop in matrix.index.values:
            rndm_dict[gene][pop] = float(0)
    for index, row in matrix.iterrows():
        mutations = int(sum(row.values))
        sample = np.random.choice(genes, p=pmf, size = mutations, replace=True)
        count_sample = Counter(sample)
        for x in count_sample.most_common():
            rndm_dict[x[0]][index] = float(x[1])

    rndm_df = pd.DataFrame.from_dict(rndm_dict, orient='index').T
    rndm_df = rndm_df.loc[:, (rndm_df != 0).any(axis=0)]
    return(rndm_df)

def calc_dispersion(Y):
    print("sad")

def get_null_results(sims = 10):
    gene_path = mydir + 'data/reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    gene_table = pd.read_csv(gene_path, sep = ' ', header = 'infer')
    gene_by_pop_path = mydir + 'data/gene_by_sample/B_S/sample_by_gene.txt'
    df = pd.read_csv(gene_by_pop_path, sep = '\t', header = 'infer', index_col = 0)
    # remove intergenic regions
    cols = [c for c in df.columns if '/' not in c]
    df = df[cols]
    df = df[['_D100' in s for s in df.index]]
    df = df[(df.T != 0).any()]
    df_rndm = get_random_gene_by_pop(df, gene_table)
    df_rndm_G = get_sample_by_gene_matrix_gscore(df_rndm)

    pcoa = get_pcoa(df_rndm)
    pcoa_xy = pcoa[0].iloc[:,0:2]
    pcoa_xy = pcoa_xy[['L0B' in s for s in pcoa_xy.index]]
    mean_x = np.mean(pcoa_xy.iloc[:,0].values)
    mean_y = np.mean(pcoa_xy.iloc[:,1].values)

    print(mean_y)

    # just L0B lines




def plot_pcoa():
    df = pd.read_csv(mydir + 'data/gene_by_sample/B_S/sample_by_gene_Gscore.txt', sep = '\t', index_col =0)
    df = df[['_D100' in s for s in df.index]]
    Y = get_pcoa(df)

    L0B = Y[['L0B' in s for s in Y.index]]
    L1B = Y[['L1B' in s for s in Y.index]]
    L2B = Y[['L2B' in s for s in Y.index]]
    L0S = Y[['L0S' in s for s in Y.index]]
    L1S = Y[['L1S' in s for s in Y.index]]
    L2S = Y[['L2S' in s for s in Y.index]]

    fig = plt.figure()
    #['#87CEEB', '#FFA500', '#FF6347']
    plt.scatter(L0B.ix[:,0].values, L0B.ix[:,1].values, marker = "o", c = '#87CEEB', s = 80)
    plt.scatter(L1B.ix[:,0].values, L1B.ix[:,1].values, marker = "o", c = '#FFA500', s = 80)
    plt.scatter(L2B.ix[:,0].values, L2B.ix[:,1].values, marker = "o", c = '#FF6347', s = 80)

    plt.scatter(L0S.ix[:,0].values, L0S.ix[:,1].values, marker = "o", s = 80, facecolors='none', edgecolors='#87CEEB')
    plt.scatter(L1S.ix[:,0].values, L1S.ix[:,1].values, marker = "o", s = 80, facecolors='none', edgecolors='#FFA500')
    plt.scatter(L2S.ix[:,0].values, L2S.ix[:,1].values, marker = "o", s = 80, facecolors='none', edgecolors='#FF6347')


    fig.tight_layout()
    fig.savefig(mydir + 'figs/test_pcoa.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

#merge_B_S_sample_by_gene_matrix()
#random_gene_by_pop()
get_null_results()
