from __future__ import division
import os, sys, copy, itertools, random, math
import numpy as np
from collections import Counter
from itertools import combinations

import scipy.stats as stats
import pandas as pd

import phylo_tools as pt


import parse_file
import timecourse_utils
import mutation_spectrum_utils
import phylo_tools as pt










#taxa = pt.taxa
taxa = ['B', 'D', 'C', 'J', 'F', 'P']
treatments = pt.treatments
replicates = pt.replicates

output_filename = pt.get_path() + '/data/fixed_gene_table.txt'
file = open(output_filename, "w")

file.write(", ".join(['Taxon', 'Treatment', 'Locus tag', 'Number fixed mutations', 'Number mutations f_max > 0.8', 'Annotation']))
file.write("\n")


for taxon in taxa:
    sys.stderr.write("Calculating gene sizes for %s...\t" % taxon)
    gene_size_map = parse_file.create_gene_size_map(taxon)
    sys.stderr.write("Done!\n")

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data


    annotation_dict =  parse_file.get_gene_annotation_dict(taxon)

    #print(gene_size_map.keys())

    # get genbank file

    for treatment in treatments:
        print(taxon, treatment)
        if treatment+taxon == '1J':
            continue
        populations = [treatment+ taxon + replicate for replicate in replicates]

        sys.stderr.write('Calculating convergence matrix for mutations...\n')

        excluded_types = set(['synonymous'])

        identifier_size_map = gene_size_map

        def get_identifier_name(gene_name):
            return gene_name

        convergence_matrix = {}
        for identifier_name in sorted(identifier_size_map.keys()):

            length = identifier_size_map[identifier_name]

            convergence_matrix[identifier_name] = {'n_muts': 0, 'n_fixed': 0}

        position_snv_dict = {}

        #all_poly_list = []
        #all_fixed_list = []
        #overlap_poly_count = 0
        #overlap_fixed_count = 0
        #polymorphic_count = 0
        #fixed_count = 0

        for population in populations:
            sys.stderr.write("Processing %s...\t" % population)

            if population in pt.populations_to_ignore:
                continue

            # calculate mutation trajectories
            # load mutations
            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)

            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple

            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            num_processed_mutations = 0

            for mutation_idx in range(0,len(mutations)):

                position, gene_name, allele, var_type, codon, position_in_codon, AAs_count, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                #Ls = haplotype_trajectories[mutation_idx]

                state_Ls = state_trajectories[mutation_idx]

                if gene_name=='intergenic':
                    continue

                if var_type in excluded_types:
                    continue

                identifier = get_identifier_name(gene_name)

                if identifier==None:
                    sys.stderr.write("No identifier for %s!\n" % gene_name)
                    continue


                good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]
                masked_state_Ls = state_Ls[good_idxs]
                #masked_Ls = Ls[good_idxs]

                if (sum([masked_state_Ls==parse_file.POLYMORPHIC][0]) == 0) and (sum([masked_state_Ls==parse_file.FIXED][0]) == 0):
                    continue



                num_processed_mutations += 1

                if sum(masked_freqs>0.8) > 0:

                    convergence_matrix[identifier]['n_muts'] += 1


                if parse_file.FIXED in masked_state_Ls:

                    convergence_matrix[identifier]['n_fixed'] += 1

            sys.stderr.write("processed %d mutations!\n" % num_processed_mutations)


        convergence_matrix_no_zeros = {k:v for k,v in convergence_matrix.items() if v['n_fixed'] > 0}


        for k, v in convergence_matrix_no_zeros.items():

            annotation_k = annotation_dict[k]
            genus = pt.genus_dict[taxon]

            if isinstance(k, list):

                locus_tag = k[0]
            else:
                locus_tag = k


            file.write(", ".join([genus, '%s-day' % str(int(10**int(treatment))), locus_tag, str(v['n_fixed']), str(v['n_muts']), annotation_k]))
            file.write("\n")



file.close()
