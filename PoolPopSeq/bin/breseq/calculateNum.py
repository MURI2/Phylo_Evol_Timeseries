from __future__ import division
import os, math, decimal, sys, argparse
#print "before scipy"
#from scipy import special
#print "after scipy"



class PPS_WT_denominator:

    '''
    This estimator of Watterson's theta was obtained from Population genomics
    from pool sequencing Ferretti et al. (doi: 10.1111/mec.12522).
    This code used here calculates the maximum composite likelihood form of
    Watterson's theta that takes into account site-specific coverage and
    corrects for ascertainment bias due the error rate.

    Here we set the error rate as 0.01, the approximate error rate for
    Illumina sequencing.
    '''

    def __init__(self, coverage, MAF = 0.01, n_c = 50):
        self.coverage = coverage
        self.MAF = MAF
        self.n_c = n_c
        # assuming population size is 10,000 (N = 10000), a convervative estimate
        # MAF = minimum allele frequency

    def harmonic_number(self, j):
        a_j = 0
        for k in range(1, j):
            a_j += 1/k
        return a_j

    def binom(self, n, k):
        if k == n:
            return 1
        elif k == 1:         # see georg's comment
            return n
        elif k > n:          # will be executed only if y != 1 and y != x
            return 0
        else:                # will be executed only if y != 1 and y != x and x <= y
            a = math.factorial(n)
            b = math.factorial(k)
            c = math.factorial(n-k)  # that appears to be useful to get the correct result
            div = a // (b * c)
            return div


    def Stirling_number_2nd(self, n, k):
        '''Stirling numbers of the second kind'''
        k_sum = 0
        normalize = decimal.Decimal(1 / math.factorial(k))
        # add up from zero to k
        for j in range(0, k+1):
            p1 = decimal.Decimal((-1) ** (k-j))
            #p2 = decimal.Decimal(special.binom(k, j))
            p2 = decimal.Decimal(self.binom(k, j))
            p3 = decimal.Decimal(j ** n)
            k_sum += p1 * p2 * p3
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

    def c_e_S(self, i):
        m_round = int(round(self.MAF * i))
        sum_1 = 0
        if m_round == 0:
            return 0
        else:
            for l in range(1, m_round + 1):
                # go up to n_c - 1
                sum_2 = 0
                for k in range(1, self.n_c):
                    p1 = (k / self.n_c) ** (l - 1)
                    p2 = (1 - (k / self.n_c)  ) ** (i - l - 1)
                    sum_2 += p1 * p2

                #sum_1 += special.binom(i, l) * sum_2
                sum_1 += self.binom(i, l) * sum_2
            return sum_1 * (1/self.n_c)

    def calculate_corrected_WT(self):
        numer_sum = 0
        for j in range(2, min(self.coverage, self.n_c) + 1):
            p1 = self.P_c(self.coverage, j) * decimal.Decimal(self.harmonic_number(j))
            p2 = decimal.Decimal(self.c_e_S(self.coverage))
            numer_sum_j = p1 - p2

            numer_sum += numer_sum_j
        return numer_sum


class PPS_pi_denominator:
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

    def __init__(self, coverage, MAF = 0.01, n_c = 50):
        self.coverage = coverage
        self.MAF = MAF
        self.n_c = n_c
        # assuming population size is 10,000 (N = 10000), a convervative estimate
        # MAF = minimum allele frequency

    def binom(self, n, k):
        if k == n:
            return 1
        elif k == 1:         # see georg's comment
            return n
        elif k > n:          # will be executed only if y != 1 and y != x
            return 0
        else:                # will be executed only if y != 1 and y != x and x <= y
            a = math.factorial(n)
            b = math.factorial(k)
            c = math.factorial(n-k)  # that appears to be useful to get the correct result
            div = a // (b * c)
            return div



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
        if m_round == 0:
            return 0
        else:
            sum_l = 0
            for l in range(1, m_round + 1):
                sum_l_k = 0
                for k in range(1, self.n_c):
                    p1 = (k / self.n_c) ** (l - 1)

                    p2 = (1 - (k / self.n_c) ) ** (i - l - 1)
                    sum_l_k += (p1 * p2)
                #sum_l += special.binom(i - 2, l - 1) * sum_l_k
                #print type(self.binom(i - 2, l - 1))
                #print self.binom(i - 2, l - 1)
                sum_l += self.binom(i - 2, l - 1) * sum_l_k
            # In SI-14 the sum of sums is multiplied by a factor of 2, which can
            # result in a negative denomenator, giving a negative result of pi
            # all the other estiators multiply the sum of sums by the inverse of the
            # number of sampled chromosomes. I need to e-mail the author for
            # clarification
            #return sum_l * 2
            return sum_l * (1/self.n_c)

    def V_i_m(self, i):
        m_round = int(float(self.MAF * i))
        sum_V_i_m_1 = 0
        for l in range(1, m_round ):
            numer = l * (i - l)
            denom = i * (i - 1)
            term_l_1 = numer / denom
            #term_l_2 = special.binom(i - 2, l - 1)
            term_l_2 = self.binom(i - 2, l - 1)
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
        E_j_m = self.E_i_m(self.coverage)
        V_j_m = self.V_i_m(self.coverage)
        return ((E_j_m ** 2) / V_j_m )



def calculateNum(mydir, OUT, n_c, MAF):
    OUT = open(OUT, 'w')
    OUT = open('/Users/WRShoemaker/Desktop/Sample_L1P1_test.txt', 'w')
    treatments = ['0', '1', '2']
    strains = ['B', 'C', 'D', 'F', 'P']
    #path = '/Users/WRShoemaker/Desktop/Sample_L1P1.txt'
    #strains = ['P']
    #treatments = ['1']
    #reps = ['1']
    reps = ['1', '2', '3', '4', '5']
    for strain in strains:
        for treatment in treatments:
            for rep in reps:
                #print strain, treatment, rep
                sample_test = 'Sample_L' + treatment + strain + rep
                path = mydir + '/' +  sample_test + '/data/coverage_merged.txt'
                if os.path.exists(path) == True:
                    #continue
                    sum_denom_pi = 0
                    sum_denom_WT = 0
                    L = 0
                    total_cov = 0
                    for line in open(path):
                        #print line
                        line = line.split()
                        L_j = int(line[2]) - int(line[1]) + 1
                        n_r_j = int(line[3])
                        # get the sites one at a time
                        if n_r_j > 2:
                            sum_denom_pi_j = PPS_pi_denominator(n_r_j, n_c = n_c, MAF = MAF).calculate_corrected_pi()
                            sum_denom_pi += sum_denom_pi_j * L_j
                            sum_denom_WT_j = PPS_WT_denominator(n_r_j, n_c = n_c, MAF = MAF).calculate_corrected_WT()
                            sum_denom_WT += sum_denom_WT_j * L_j
                        L += L_j
                        total_cov += n_r_j
                    mean_cov = total_cov / L
                    print>> OUT, sample_test, L, mean_cov, sum_denom_pi, sum_denom_WT

    OUT.close()

#if __name__ == '__main__':
#    parser = argparse.ArgumentParser(description = "Calculate denominator for pop gen estimators")
#    parser.add_argument('-n', type = int, default = 70, help = "number of chromosomes")
#    parser.add_argument('-m', type = float, default = 0.01, help = "minimum acceptable allele frequency")
    #parser.add_argument('-o', type = str, default = "", help = "out file")
    #OUT = "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/D100/sample_numerators.txt"
OUT = "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/D100/D100_denominator.txt"
mydir = os.path.expanduser("/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/D100")

#test_in = ''
#test_out = ''

#params = parser.parse_args()
#n_c = params.n
#MAF = params.m

n_c = 70
MAF = 0.01
calculateNum(mydir, OUT, n_c, MAF)

#calculateNum(mydir, OUT, n_c, MAF)
