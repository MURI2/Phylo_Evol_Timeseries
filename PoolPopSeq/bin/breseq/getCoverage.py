from __future__ import division
import os, math, decimal, sys, argparse
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def condenseCoverage(IN, OUT_path):
    OUT = open(OUT_path, 'w')
    with open(IN) as f:
        coverage_max = 0
        position_start = 0
        count = 0
        for x ,line in enumerate(f):

            line_list = line.split()
            coverage = line_list[3]
            contig = line_list[0]
            position = line_list[1]
            if x == 0:
                coverage_max = coverage
                position_start = position
                count == 0
            else:
                if coverage == coverage_max:
                    count += 1
                else:
                    print>> OUT, contig, position_start, int(position_start) + int(count), coverage_max

                    coverage_max = coverage
                    position_start = position
                    count = 0

        print>> OUT, contig, position_start, int(position_start) + int(count), coverage_max

    OUT.close()

def mean_std_window(IN, OUT_path, block_size = 100):
    OUT = open(OUT_path, 'w')
    print>> OUT, 'contig', 'start', 'stop', 'cov_mean', 'cov_std'
    with open(IN) as f:
        coverage_list = []
        position_list = []
        contig_list = []
        for x ,line in enumerate(f):
            line_list = line.split()
            coverage = line_list[3]
            contig = line_list[0]
            position = line_list[1]
            coverage_list.append(int(coverage))
            position_list.append(int(position))
            contig_list.append(contig)
            if  x % block_size == 0 and x != 0:
                contig_set = set(contig_list)
                if len(contig_set) > 1:
                    contig_out = ','.join(list(contig_set))
                else:
                    contig_out = list(contig_set)[0]
                print>> OUT, contig_out, position_list[0], position_list[-1], \
                        str(np.mean(coverage_list)), str(np.std(coverage_list))
                coverage_list = []
                position_list = []
                contig_list = []
    OUT.close()

def mean_std_window_plot(OUT_path, OUT_plot_path):
    cov_df = pd.read_csv(OUT_path, sep = ' ')
    x = cov_df.start.values
    y = cov_df.cov_mean.values
    y_erros = cov_df.cov_std.values * 2
    fig = plt.figure()
    ax = fig.add_subplot(1, 1,  1)
    ax.plot(x, y, lw = 2, color = '#FF6347')
    ax.fill_between(x, y+y_erros, y-y_erros, facecolor='#FF6347', alpha=0.5)
    plt.axhline(y = 100, c = 'grey', linestyle = '--', lw = 3)
    ax.set_xlabel('Posistion (100 bp)', fontsize=20)
    ax.set_ylabel('Coverage', fontsize=14)
    ax.set_ylim([0, 120])
    ax.set_xlim([x[0], x[-1]])
    fig.tight_layout()
    fig.savefig(OUT_plot_path + '.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def selectSites(OUT1_path, OUT2_path):
    poly_dict = {}
    #OUT_subset =
    OUT_groups = OUT2_path.split('/')
    OUT_subset = open('/'.join(OUT_groups[:-1]) + '/coverage_merged_subset.txt', 'w')
    with open(OUT2_path) as f:
        for x ,line in enumerate(f):
            line_split = line.split()
            if len(line_split) < 3:
                continue
            seq_id = line_split[3]
            position = line_split[4]
            if seq_id in poly_dict:
                poly_dict[seq_id].append(int(position))
            else:
                poly_dict[seq_id] = [int(position)]
    with open(OUT1_path) as g:
        # g = gd file
        contig = ''
        L = 0
        for x ,line in enumerate(g):
            if line_split[0] != contig:
                contig = line_split[0]
            line_split = line.split()

            contig = line_split[0]
            sites = poly_dict[contig]
            l2 = [i for i in sites if i >= int(line_split[1]) and i <= int(line_split[2])]
            if len(l2) > 0:
                for l in l2:
                    sites.remove(l)
                    print int(line_split[1]), l, int(line_split[2])
                    print>> OUT_subset, L, contig, line_split[1], line_split[2], line_split[3]
                    print len(sites)
                L += 1

    OUT_subset.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "condense coverage")
    parser.add_argument('-c', action='store_true', default=False)
    parser.add_argument('-m', action='store_true', default=False)
    parser.add_argument('-v', action='store_true', default=False)

    parser.add_argument('-i', type = str, default = "", help = "in file")
    parser.add_argument('-o', type = str, default = "", help = "out file")
    #parser.add_argument('-p', type = str, default = "", help = "out plot")

    params = parser.parse_args()

    if params.v == False and params.c == True:
        condenseCoverage(params.i, params.o)
    elif params.v == False and params.m == True:
        mean_std_window(IN = params.i, OUT_path = params.o)
    #    if params.p != "":
    #        mean_std_window_plot(OUT_path = params.o, OUT_plot_path = params.p)
    else:
        print "No argument provided\nExiting program"
