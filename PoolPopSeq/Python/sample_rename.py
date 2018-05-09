from __future__ import division
import os, re
from subprocess import call

runs = ['161213', '170303', '170623', '170721']

def rename():
    for run in runs:
        for i, x in enumerate(list(os.walk('/N/dc2/projects/muri2/Task2/PoolPopSeq/data/' + run))):
            if i  == 0:
                continue
            path = x[0]
            for y in x[2]:
                if y == 'SampleSheet.csv':
                    continue
                if run == '161213':
                    y_split = y.split('_', 1)
                    y_out = y_split[0] + '-100_' + y_split[1].split('.', 1)[0] + '_' + \
                        run + '.' + y_split[1].split('.', 1)[1]
                elif run == '170303':
                    y_out = y.split('.', 1)[0] + '_' + run + '.' + y.split('.', 1)[1]
                elif run == '170721':
                    y_rep = y.replace('-D10', '')
                    y_rep2 = y_rep.replace('OF3', '0F3')
                    y_out = 'L' + y_rep2.split('.', 1)[0] + '_' + run + '.' + y_rep2.split('.', 1)[1]
                elif run == '170623':
                    y_split = y.split('_', 4)
                    y_out = 'L' + y_split[1] + y_split[2] + '-' + y_split[3][1:] \
                        + '_' + y_split[4].split('.', 1)[0] + '_' + run + '.' + y_split[4].split('.', 1)[1]
                y_out_path1 = '/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/' \
                        + 'D' +  re.split('-|_',y_out)[1] + '/'
                if os.path.exists(y_out_path1) != True:
                    os.makedirs(y_out_path1)
                y_out_path2 = y_out_path1 + 'Sample_' + y_out.split('_', 1)[0] + '/'
                if os.path.exists(y_out_path2) != True:
                    os.makedirs(y_out_path2)
                y_out_path = y_out_path2 + y_out
                y_path = x[0] + '/' + y
                call(["cp", y_path, y_out_path])
                print y_path
                print y_out_path


#L0B2-100_GTAGAGGA-AAGGAGTA_L002_R1_001.fastq.gz
#L2P2-300_TCGACGTC-TATCCTCT_L002_R1_001_170623.fastq.gz

rename()
