from __future__ import division
import os, math, numbers, itertools, re
import pandas as pd
import numpy as np
from string import maketrans
from collections import Counter

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")



class cleanBreseq_html:

    def __init__(self, path):
        self.path = path
        self.IN = pd.read_html(self.path)
        self.freqs = self.IN[1]
        self.freqs.columns = self.freqs.iloc[1]
        self.freqs = self.freqs.drop(self.freqs.index[[0,1]])

    def percent_to_freq(self, x):
        percent = float(x.replace('%', ''))
        return percent / 100

    def arrowsToASCII(self, x, column):
        minus = u"\u2011"
        R_arrow = u"\u2192"
        L_arrow = u"\u2190"
        delta = u"\u0394"
        if column == 'gene':
            x_split = x.split()
            if len(x_split) == 1:
                out = x_split[0]
            elif len(x_split) == 2:
                if x_split[1] == u"\u2192":
                    arrow = '->'
                elif x_split[1] == u"\u2190":
                    arrow = '<-'
                else:
                    pass
                gene_name = x_split[0].replace(minus, '-')
                out = gene_name + '_' + arrow
            else:
                if x_split[1] == u"\u2192":
                    L_arrow = '->'
                elif x_split[1] == u"\u2190":
                    L_arrow = '<-'
                if len(x_split) > 3:
                    if x_split[3] == u"\u2192":
                        R_arrow = '->'
                    elif x_split[3] == u"\u2190":
                        R_arrow = '<-'
                    elif x_split[3] == u"\u2013":
                        R_arrow = '-'
                gene_name1 = x_split[0].replace(minus, '-')
                if len(x_split) > 4:
                    gene_name2 = x_split[4].replace(minus, '-')
                    out =  gene_name1 + '_' + L_arrow + '_' + R_arrow + '_' + gene_name2
                elif len(x_split) > 3:
                    out =  gene_name1 + '_' + L_arrow + '_' + R_arrow
                else:
                    out =  gene_name1 + '_' + L_arrow
            out = ''.join([i if ord(i) < 128 else '' for i in out])
            return out
        elif column == 'annotation':

            if isinstance(x, numbers.Number) == True:
                out = ''
                return out
            else:
                x_split = x.split()
                if len(x_split) == 2:
                    x_1_split = x_split[1].split()
                    if minus in x_1_split[0]:
                        out = x_1_split[0].replace(minus, '-')
                        out = out.replace('(', '')
                        out = out.replace(')', '')
                        #out = out.replace('/', '_')
                    elif R_arrow in x_1_split[0]:
                        out = x_1_split[0].replace(R_arrow, '->')
                        out = out.replace('(', '')
                        out = out.replace(')', '')
                    else:
                        out = x_1_split[0]
                    out = out.replace(u"\u2013", '')
                    #out = ''.join([i if ord(i) < 128 else '' for i in out])
                    return out
                elif len(x_split) == 3:
                    x_1_split = x_split[1].split()

                    middle = x_split[1].replace(minus, '_')
                    middle = middle.replace('(', '_')
                    end = x_split[2].replace(')', '')
                    out = x_split[0] + middle + '_' + end
                    out = ''.join([i if ord(i) < 128 else '' for i in out])
                    return out
                else:
                    pass

        elif column == 'mutation':
            x_split = x.split()
            if len(x_split) == 1:
                out = x_split[0].replace(R_arrow, '->')
                out = out.replace(L_arrow, '<-')
                out = ''.join([i if ord(i) < 128 else '' for i in out])
            elif len(x_split) == 2:
                out1 = x_split[0].replace(delta, 'delta_')
                out2 = x_split[1].replace(R_arrow, '->')
                out2 = out2.replace('bp', '_bp')
                out = out1 + out2
                out = out.replace(',', '')
                out = ''.join([i if ord(i) < 128 else '' for i in out])
            else:
                out = '_'.join(x_split)
                out = out.replace(R_arrow, '->')
            return out

        elif column == 'description':
            out = x.replace(u"\u2011", '-')
            out = x.replace(u"\u2013", '-')
            out = x.replace(' ', '_')
            out = x.replace(' ', '_')
            #out = x.replace(',', '')
            out = out.replace(minus, '-')
            out = ''.join([i if ord(i) < 128 else '' for i in out])

            return out
        elif column == 'evidence':
            to_remove = '\xa0'
            #out = x.replace(to_remove, '_')
            #  #out = ''.join([i if ord(i) < 128 else '' for i in out])
            x_split = x.split()
            if len(x_split) == 1:
                out = x_split[0]
            else:
                out = x_split[0] + '_' + x_split[0]
            out = ''.join([i if ord(i) < 128 else '' for i in out])
            return out

    def getCoverage(self, evidence_path):
        coverage_list = []
        for index, row in self.freqs.iterrows():
            mutation = row['mutation']
            mutation_number = str(index - 1)
            if '->' in mutation and 'bp' not in mutation:
                prefix = 'SNP_'
                path = evidence_path + prefix + mutation_number + '.html'
                if '/D100/Sample_L0C4/output/evidence/SNP_162.html' in path:
                    coverage = float('NaN')
                elif os.path.exists(path) == True:
                    IN = pd.read_html(path)
                    IN_1 = IN[1]
                    coverage = IN_1.iloc[2][8]
                else:
                    coverage = float('nan')
            elif '->' in mutation and 'bp' in mutation :
                prefix = 'SUB_'
                IN = pd.read_html(evidence_path + prefix + mutation_number + '.html')
                IN_1 = IN[1]
                coverage = IN_1.iloc[2][8]
                coverage_split = coverage.split()
                if len(coverage_split) > 1:
                    coverage = coverage_split[1]
                else:
                    coverage = coverage
            elif 'delta' in mutation:
                prefix = 'DEL_'
                path = evidence_path + prefix + mutation_number + '.html'
                if os.path.exists(path) == True:
                    IN = pd.read_html(path)
                    IN_1 = IN[1]
                    coverage = IN_1.iloc[2][8]
                else:
                    coverage = float('nan')
            elif '+' in mutation:
                prefix = 'INS_'
                path = evidence_path+ prefix + mutation_number + '.html'
                if os.path.exists(path) == True:
                    IN = pd.read_html(path)
                    IN_1 = IN[1]
                    if '%' in IN_1.iloc[2][8]:
                        coverage = float('nan')
                    else:
                        coverage = IN_1.iloc[2][8]
                else:
                    coverage = float('nan')
            else:
                coverage = float('nan')

            try:
                float(coverage)
            except ValueError:
                if '%' in coverage:
                    coverage = float('nan')
            coverage_list.append(coverage)

        return coverage_list


    def cleanColumns(self, evidence_path):
        self.freqs['gene'] = self.freqs['gene'].apply(lambda x: self.arrowsToASCII(x, 'gene'))
        self.freqs['annotation'] = self.freqs['annotation'].apply(lambda x: self.arrowsToASCII(x, 'annotation'))
        self.freqs['mutation'] = self.freqs['mutation'].apply(lambda x: self.arrowsToASCII(x, 'mutation'))
        self.freqs['description'] = self.freqs['description'].apply(lambda x: self.arrowsToASCII(x, 'description'))
        self.freqs['evidence'] = self.freqs['evidence'].apply(lambda x: self.arrowsToASCII(x, 'evidence'))
        self.freqs['freq'] = self.freqs['freq'].apply(lambda x: self.percent_to_freq(x))

        if self.freqs.shape[1] == 8:
            self.freqs.columns = ['evidence', 'seq_id', 'position', 'mutation', 'freq', 'annotation', 'gene', 'description']
            self.freqs['seq_id_position'] = self.freqs['seq_id'].astype(str) + '_' +self.freqs['position'].astype(str)
        coverage = pd.Series(self.getCoverage(evidence_path))
        self.freqs['coverage'] = coverage.values
        return self.freqs



def cleanSamples_html():
    #treatments = ['0', '1', '2']
    strains = ['B', 'C', 'D', 'F', 'J', 'P']
    #strains = ['F', 'J', 'P']
    #strains = ['P']
    treatments = ['1', '2']
    reps = ['1', '2', '3', '4', '5']
    #reps = ['1', '2', '3', '4', '5']
    for treatment in treatments:
        for strain in strains:
            for rep in reps:
                path = mydir+ 'breseq_output_gbk/D100/Sample_L' + treatment + strain + rep + '/output/index.html'
                #path = mydir+ 'output/index.html'
                if os.path.exists( path) == True:
                    print treatment + strain + rep
                    evidence_path = mydir + 'breseq_output_gbk/D100/Sample_L' + treatment + strain + rep + '/output/evidence/'
                    clean = cleanBreseq_html(path).cleanColumns(evidence_path)
                    OUTname = mydir + 'breseq_output_gbk_clean/D100/Sample_L' + treatment + strain + rep + '.txt'
                    clean.to_csv(OUTname, sep = '\t', index = False)


def uniqueMutations_html():
    #treatments = ['0', '1', '2']
    #strains = ['B', 'C', 'D', 'F', 'J', 'P']
    strains = ['D']
    #reps = ['1', '2', '3', '4', '5']
    treatments = ['0']
    reps = ['1', '2']
    for strain in strains:
        count = 0
        for treatment in treatments:
            for rep in reps:
                path =  mydir + 'breseq_output_gbk_clean/D100/Sample_L' +  treatment + strain + rep + '.txt'
                if os.path.exists(path) == True:
                    IN = pd.read_csv(path, sep = '\t', header = 'infer')
                    freq = 'freq_L' + treatment +  strain + rep
                    IN = IN.rename(columns = {'freq':freq})
                    coverage = 'coverage_L' + treatment +  strain + rep
                    IN = IN.rename(columns = {'coverage':coverage})
                    if count == 0:
                        merged = IN
                        merged_freq = merged[freq]
                        merged.drop(labels=[freq], axis=1,inplace = True)
                        merged.insert(len(merged.columns)-1, freq, merged_freq)
                    else:
                        if 'seq_id' in IN.columns and 'seq_id_position' in IN.columns:
                            merged = pd.merge(merged, IN[[freq, coverage, 'position', 'seq_id', 'seq_id_position']], \
                                how='outer', on=['position', 'seq_id', 'seq_id_position'])
                        else:
                            merged = pd.merge(merged, IN[[freq, coverage, 'position']], \
                                how='outer', on=['position'])
                    count += 1
        OUTname = mydir + 'breseq_output_gbk_clean_merged/D100/Strain_' + strain + '.txt'
        if 'seq_id' in merged.columns and 'seq_id_position' in merged.columns:
            merged = merged.drop_duplicates(subset = ['position', 'seq_id', 'seq_id_position'])
        else:
            merged = merged.drop_duplicates(subset = ['position'])
        merged.to_csv(OUTname, sep = '\t', index = False)

        sample_freqs = [x for x in merged.columns if 'freq_' in x]
        merged_NoDups = merged[merged[sample_freqs].apply(lambda x: min(x) != max(x), 1)]
        OUTnameNoDups = mydir + 'breseq_output_gbk_clean_merged_NoDups/D100/Strain_' + strain + '.txt'
        merged_NoDups.to_csv(OUTnameNoDups, sep = '\t', index = False)
