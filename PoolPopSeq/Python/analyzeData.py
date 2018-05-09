from __future__ import division
import os, math, decimal
import pandas as pd
import numpy as np
from scipy import special

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")

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

species_dict = {'B': 'Bacillus', 'C':'Caulobacter', 'D':'Deinococcus', \
    'F': 'Pedobacter', 'J': 'Janthinobacterium', 'P':'Pseudomonas', \
    'S': 'Bacillus_spoA'}


class PPS_Wtheta:

    '''
    This estimator of Watterson's theta was obtained from Population genomics
    from pool sequencing Ferretti et al. (doi: 10.1111/mec.12522).
    This code used here calculates the maximum composite likelihood form of
    Watterson's theta that takes into account site-specific coverage and
    corrects for ascertainment bias due the error rate.

    Here we set the error rate as 0.01, the approximate error rate for
    Illumina sequencing.
    '''

    def __init__(self, freqs, coverage, MAF = 0.01, n_c = 70):
        self.freqs = freqs
        self.coverage = coverage
        self.MAF = MAF
        self.n_c = n_c
        # assuming population size is 10,000 (N = 10000), a convervative estimate
        # MAF = minimum allele frequency

    def remove_MAFs(self):
        filtered_freqs = [x for x in self.freqs if x >= self.MAF]
        return filtered_freqs

    def harmonic_number(self, j):
        a_j = 0
        for k in range(1, j):
            a_j += 1/k
        return a_j

    def Stirling_number_2nd(self, n, k):
        '''Stirling numbers of the second kind'''
        if n == k:
            return 1
        k_sum = 0
        normalize = decimal.Decimal(1 / math.factorial(k))
        #normalize = 1 / math.factorial(k)
        # add up from zero to k
        for j in range(0, k+1):
            S_p1 = decimal.Decimal((-1) ** (k - j))
            S_p2 = decimal.Decimal(special.binom(k, j))
            S_p3 = decimal.Decimal(j ** n)
            #S_p1 = int((-1) ** (k - j))

            #S_p2 = int(special.binom(k, j))
            #S_p3 = int(j ** n)
            k_sum += S_p1 * S_p2 * S_p3

        return k_sum * normalize


    def P_c(self, n_r, j):
        '''
        Calculates the probability that a set of n_r sequences randomly
        extracted (with repetitions) from n_c possible lineages contains
        sequences coming from precisely j different lineages

        n_r = read depth at a site
        n_c = number of independant lineages (i.e. chromosomes) in the pool
        '''
        stirling = self.Stirling_number_2nd(n_r, j)
        numer = decimal.Decimal(stirling) * decimal.Decimal(math.factorial(self.n_c))
        denom = decimal.Decimal(math.factorial(self.n_c -j)) * decimal.Decimal(self.n_c ** n_r)
        return numer / denom

    def c_e_S(self, n_r):
        m_round = int(round(self.MAF * n_r))
        if m_round == 0:
            return 0
        else:
            sum_l = 0
            for l in range(1, m_round + 1):
                # go up to n_c - 1
                sum_l_k = 0
                for k in range(1, self.n_c):
                    p1 = (k / self.n_c) ** (l - 1)
                    p2 = (1 - (k / self.n_c)  ) ** (n_r - l - 1)
                    sum_l_k += p1 * p2

                sum_l += special.binom(n_r, l) * sum_l_k
            return sum_l * (1/self.n_c)

    def calculate_corrected_WTheta(self):
        freqs = self.remove_MAFs()
        S = len(freqs)
        L = 0
        # you should be iterating over ALL sites, not just sites with a variant.
        # this takes into acocunt the length of the genome (similar to dividing by L)
        numer_sum = 0
        for i in self.coverage:
            numer_sum_i = 0
            L_i = i[1] - i[0] + 1
            # reads (i.e. coverage)
            n_r_i = i[2]
            for j in range(2, min(n_r_i, self.n_c) + 1):
                P_c = self.P_c(n_r_i, j)
                a_j = decimal.Decimal(self.harmonic_number(j))
                #c_e_S = decimal.Decimal(self.c_e_S(n_r_i))
                #p2 = 0
                numer_sum_i += P_c * a_j
            c_e_S = decimal.Decimal(self.c_e_S(n_r_i))
            numer_sum_i = numer_sum_i - c_e_S

            numer_sum += numer_sum_i * L_i
            L += L_i
        #print S /  (self.harmonic_number(self.n_c)  * L)
        return S / numer_sum

class PPS_pi:
    '''
    This estimator of Tajima's theta/ nucleotide diversity/
    was obtained from Population genomics from pool sequencing
    Ferretti et al. (doi: 10.1111/mec.12522).
    This code used here calculates the maximum composite likelihood form of
    pi that takes into account site-specific coverage and
    corrects for ascertainment bias due the error rate.

    Here we set the error rate as 0.01, the approximate error rate for
    Illumina sequencing.

    Because this is pooled sequencing pi is calculated as the site-specific
    diversity for a bi-allelic site. See below URL for an example
    got here: http://genapps.uchicago.edu/slider/index.html
    Also in Walsh & Lynch 201?
    '''

    def __init__(self, freqs, coverage, MAF = 0.01, n_c = 70, better_estimate = True):
        self.freqs = freqs
        self.coverage = coverage
        self.MAF = MAF
        self.n_c = n_c
        self.better_estimate = better_estimate
        self.freqs_cov = self.remove_MAFs()
        # assuming population size is 10,000 (N = 10000), a convervative estimate
        # MAF = minimum allele frequency
        # g = genome size

    def remove_MAFs(self):
        filtered_freqs = [x for x in self.freqs if x >= self.MAF]
        return filtered_freqs

    def calculate_pi(self):
        sum_pi = 0
        for x in self.freqs_cov:
            sum_pi += ( 2 * x[0] * (1 - x[0]))
        return sum_pi

    def V_i(self, i):
        term_1_1 = (i - 2) / i
        term_1_2 = (i - 3) / (i - 1)
        term_1_3 = (self.n_c - 1) / self.n_c
        term_1_4 = (self.n_c - 2) / self.n_c
        term_1_5 = (self.n_c - 3) / self.n_c

        term_1 = (1/3) * term_1_1 * term_1_2 * term_1_3 * term_1_4 * term_1_5

        term_2_1 = (i - 2) / i
        term_2_2 = (i - 3) / (i - 1)
        term_2_3 = ((self.n_c - 1) / self.n_c ) ** 2

        term_2 = (2 / self.n_c) * term_2_1 * term_2_2 * term_2_3

        term_3 = (2 / i) * ( (self.n_c - 1) / self.n_c )
        return term_1 + term_2 + term_3

    def c_e_pi(self, i):
        # i = coverage
        # m_round = lower bound on abundance of an allele
        m_round = int(round(self.MAF * i))
        sum_l = 0
        for l in range(1, m_round + 1):
            sum_l_k = 0
            for k in range(1, self.n_c):
                p1 = (k / self.n_c) ** (l - 1)

                p2 = (1 - (k / self.n_c) ) ** (i - l - 1)
                sum_l_k += (p1 * p2)
            sum_l += special.binom(i - 2, l - 1) * sum_l_k
        # In SI-14 the sum of sums is multiplied by a factor of 2, which can
        # result in a negative denomenator, giving a negative result of pi
        # all the other estiators multiply the sum of sums by the inverse of the
        # number of sampled chromosomes. I need to e-mail the author for
        # clarification
        return sum_l * (1/self.n_c)

    def V_i_m(self, i):
        m_round = int(float(self.MAF * i))
        sum_V_i_m_1 = 0
        for l in range(1, m_round ):
            numer = l * (i - l)
            denom = i * (i - 1)
            term_l_1 = numer / denom
            term_l_2 = special.binom(i - 2, l - 1)
            sum_V_i_m_2 = 0
            for k in range(1, self.n_c):
                term_k_1 = (k / self.n_c) ** (l - 1)
                term_k_2 = (1 - (k / self.n_c)  ) ** (i - l - 1)
                sum_V_i_m_2 += term_k_1 * term_k_2
            sum_V_i_m_1 += sum_V_i_m_2 * term_l_1 * term_l_2
        return self.V_i(i) - (sum_V_i_m_1 * 4)

    def E_i_m(self, coverage):
        self.c_e_pi(coverage)
        return ((self.n_c - 1) / self.n_c) - self.c_e_pi(coverage)


    def calculate_corrected_pi(self):
        sum_pi_numer = 0
        sum_pi_denom = 0
        pi_sum_test = 0
        L = 0
        for i in self.freqs:
            # second element is coverage
            E_i_m = self.E_i_m(i[1])
            V_i_m = self.V_i_m(i[1])
            pi_i_m = 2 * i[0] * (1 - i[0])
            pi_sum_test += 2 * i[0] * (1 - i[0])
            sum_pi_numer += (pi_i_m * E_i_m) / V_i_m
        for j in self.coverage:
            # third element is coverage
            E_j_m = self.E_i_m(j[2])
            V_j_m = self.V_i_m(j[2])
            L_j = j[1] - j[0] + 1
            L += L_j
            sum_pi_denom += ((E_j_m ** 2) / V_j_m ) * L_j

        #print (pi_sum_test / L) *  (  self.n_c  / (self.n_c - 1)  )
        return sum_pi_numer / sum_pi_denom


class PPS_tajimas_D:
    '''
    A class to estimate Tajima's D using pi (Tajima's theta), S
    (number polymorphic sites), and N (pop. size)
    We'll use mean coverage here for N. I need to get a better estimate for
    the variance term of Tajima's D
    '''

    def __init__(self, pi, w, S, n):
        self.pi = float(pi)
        self.w = float(w)
        self.S = int(S)
        self.n = int(n)

    def a1(self):
        '''Given n, this function returns the (n-1)th harmonic number'''
        return sum((1.0/d) for d in range(1,self.n))

    def a2(self):
        '''Given n, this function returns the (n-1)th squared harmonic number'''
        return sum((1.0/(d**2)) for d in range(1,self.n))

    def b1(self):
        '''Creates b1 for the variance of Tajima's theta'''
        return ((self.n+1) /  (3*(self.n-1)))

    def b2(self):
        '''Creates b2 for the variance of Tajima's theta'''
        num = ((self.n**2) + self.n + 3) * 2
        den = 9 * self.n * (self.n-1)
        return num / den

    def c1(self):
        '''Creates c1 for the variance of Tajima's theta'''
        return self.b1() - (1 / self.a1())

    def c2(self):
        '''Creates c2 for the variance of Tajima's theta'''
        return self.b2() - ((self.n+2) / (self.a1() * self.n)) + (self.a2() / (self.a1() ** 2 ))

    def e1(self):
        return self.c1() / self.a1()

    def e2(self):
        return self.c2() / ( (self.a1() ** 2) + self.a2() )

    def tajimas_D(self):
        num = self.pi - self.w
        den = math.sqrt( (self.e1() * self.S) + (self.e2() * self.S * (self.S-1)) )
        T_D = num / den
        return T_D


def harmonic_number(j):
    a_j = 0
    for k in range(1, j):
        a_j += 1/k
    return a_j

def popGenTable(MAF = 0.01, n_c = 70):
    L_samples = {'C': 4042929, 'D': 3284156, 'F': 5836693, 'P':6592875, \
        'B':4215606, 'J':6082545, 'S':4215606}
    strains = ['B', 'C', 'D', 'F', 'J', 'P', 'S']
    #strains = ['B']
    OUT =  mydir + 'pop_gen_stats/D100/popGenTable.txt'
    OUT = open(OUT, 'w')
    print>> OUT, 'strain', 'treatment', 'replicate', 'k', 'S', 'pi', 'W_T', \
            'T_D', 'k_L', 'S_L', 'pi_L', 'W_T_L', 'mean_coverage'
    for strain in strains:
        path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/D100/Strain_' + strain + '_D100_SNP.txt'

        if os.path.exists(path) == True:
            IN = pd.read_csv(path, sep = '\t', header = 'infer')
            sample_freqs = [x for x in IN.columns if 'frequency_' in x]
            number_samples = len(sample_freqs)
            for sample_freq in sample_freqs:
                sample_freq_split = sample_freq.split('_')
                total_cov = 'total_cov_' + sample_freq_split[1] + '_D100'
                sample_freq_cov = list(zip( IN[sample_freq], IN[total_cov]))
                sample_freq_cov = [x for x in sample_freq_cov if np.isnan(x[0]) == False]
                k = len([x for x in sample_freq_cov if x[0] == float(1)])
                poly = [x for x in sample_freq_cov if x[0] != float(1) and x[0] > MAF]
                S = len(poly)
                W_T = S /harmonic_number(n_c)
                pi = sum([(2 * x[0] * (1-x[0])) for x in poly]) * ((n_c - 1) /  n_c)
                coverage = 0
                for site in poly:
                    site_split = site[1].split('/')
                    coverage += (int(site_split[0]) + int(site_split[1]))
                if coverage > 0:
                    mean_coverage = coverage / S
                else:
                    #print sample_freq
                    mean_coverage = 0
                if mean_coverage == 0:
                    T_D = float('nan')
                else:
                    T_D = PPS_tajimas_D(pi, W_T, S, mean_coverage).tajimas_D()
                L = L_samples[strain]
                k_L = k /L
                S_L = S / L
                pi_L = pi / L
                W_T_L = W_T / L
                treatment = sample_freq_split[1][1]
                replicate = sample_freq_split[1][3]
                print>> OUT, species_dict[strain], treatment, replicate, k, S, pi, \
                        W_T, T_D, k_L, S_L, pi_L, W_T_L, mean_coverage
    OUT.close()

def get_column_name(row, row_name):
    names = row.index.tolist()
    frequency_samples = [y for y in names if 'frequency_' in y]
    frequency_row = row[frequency_samples]
    sample = frequency_row[frequency_row.notnull()].index.tolist()[0]
    if row_name == 'treatment':
        return sample[-3]
    elif row_name == 'replicate':
        return sample[-1]
    else:
        return float('NaN')



def geneTable(MAF = 0.01):
    strains = ['B', 'C', 'D', 'F', 'J', 'P']
    #strains = ['B']
    to_rename = ['frequency', 'total_cov', 'number', 'file_number', \
        'prediction', 'consensus_score', 'polymorphism_score', \
        'fisher_strand_p_value', 'ks_quality_p_value', 'bias_e_value', \
        'bias_p_value', 'reject', 'snp_type', 'type', \
        'major_base', 'minor_base', 'sample']
    for strain in strains:
        path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/D100/Strain_' + strain + '.txt'
        if os.path.exists(path) == True:
            IN = pd.read_csv(path, sep = '\t', header = 'infer')
            gene_name = IN['gene_name']
            genes = set(gene_name.values)
            indexes = []
            for gene in genes:
                rows = IN.loc[IN['gene_name'] == gene]
                row_count = rows.shape[0]
                if row_count > 1:
                    gene_indexes = rows.index.values
                    indexes.extend(gene_indexes)

            path = mydir + 'multiple_gene_hits/D100/Strain_' + strain + '.txt'
            # find non-nan values and keep those columns. set new column with
            # strain, treatment, and replicate
            out_df = IN.iloc[indexes]
            out_df['treatment'] = out_df.apply(get_column_name, args=('treatment',), axis=1)
            out_df['replicate'] = out_df.apply(get_column_name, args=('replicate',), axis=1)
            for x in to_rename:
                x_cols = [y for y in out_df.columns if (x + '_') in y]
                # merge the columns
                out_df[x] = out_df[x_cols].fillna('').sum(axis=1)
                # remove the columsn you just merged
                out_df = out_df.drop(x_cols, axis=1)

            out_df.to_csv(path, sep = '\t', index = False)

            out_df_small = out_df[['seq_id', 'position', 'gene_name', \
                'gene_position', 'treatment', 'replicate', 'gene_product']]
            out_df_small.to_csv(mydir + 'multiple_gene_hits/D100/Strain_' \
                + strain + '_reduced.txt', sep = '\t', index = False)






#freqs = [(0.5, 100), (0.5, 100 ), (0.5, 100)]
#coverage = [80, 80, 80, 80, 80, 80, 80, 80, 80, 80]
# L = coverage[i][1] - coverage[i][0] + 1
#coverage = [(1,5,100), (6,10,100), (11,20,100)]
#print PPS_Wtheta(freqs, coverage).calculate_corrected_WTheta()
#print PPS_pi(freqs, coverage).calculate_corrected_pi()

popGenTable()
#geneTable()
