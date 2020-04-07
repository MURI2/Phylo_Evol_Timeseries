from __future__ import division
import sys, os, bz2
import parse_file
import phylo_tools as pt
import numpy as np
import timecourse_utils
import statsmodels.stats.multitest as multi
#mydir = os.path.expanduser("~/GitHub/LTEE-metagenomic-master/")

taxon = sys.argv[1]

treatments = pt.treatments
replicates = pt.replicates

gene_data = parse_file.parse_gene_list(taxon)

position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = parse_file.create_annotation_map(taxon=taxon, gene_data=gene_data)

#FDR = 0.05
FDR = 0.1
pvalue_idx = 6

total_num_passed = 0


sys.stdout.write("Filtering trajectories with FDR=%g\n" % FDR)

populations = pt.get_populations(taxon)

# first determine pvalue threshold
# across all nonmutator populations

#for population in ['0B1']:
#for treatment in ['0']:
    #pvalues = []
    #for replicate in ['1','2','3','4','5']:
    #    population = treatment + taxon + replicate
    #    #file = open(input_filename_template % population,"r")
    #    likelihood_filename = '%s_likelihood_timecourse.bz' % (population)
    #    likelihood_timecourse_path = pt.get_path() + '/data/timecourse_likelihood/' + likelihood_filename
    #    file = bz2.open(likelihood_timecourse_path,"rt")
    #    file.readline() # depth line!
    #    for line in file:

    #        items = line.split(",")
    #        location = int(items[1])
    #        pitems = items[pvalue_idx].split()
    #        pvalue = float(pitems[1])
    #        autocorrelation = float(pitems[0])

    #        if location==0:
    #            print(line)

    #        if parse_file.annotate_gene(location, position_gene_map)=='repeat':
    #            continue
    #
    #        pvalues.append(pvalue)
    #    file.close()

    #ntot = len(pvalues)

    #pvalues = np.array(pvalues)

    #for p in sorted(pvalues,reverse=True):
    #    print(p, (ntot*p/(pvalues <= p).sum()) )
    #    if ntot*p/(pvalues <= p).sum() < FDR:
    #        break

    #reject, pvals_corrected, alphacSidak, alphacBonf = multi.multipletests(pvalues, alpha=0.8, method = 'fdr_bh', is_sorted=False)
    #
    #print(alphacSidak, alphacBonf)
    # this is now the global threshold pvalue
    #nonmutator_threshold_pvalue = p
    #sys.stdout.write(treatment + taxon +" p(FDR) = %g\n" % nonmutator_threshold_pvalue)


for treatment in treatments:
    for replicate in replicates:
        population = treatment + taxon + replicate
        if population in pt.populations_to_ignore:
            continue
        #file = open(input_filename_template % population,"r")
        likelihood_filename = '%s_likelihood_timecourse.bz' % (population)
        likelihood_timecourse_path = pt.get_path() + '/data/timecourse_likelihood/' + likelihood_filename
        if os.path.exists(likelihood_timecourse_path)  == False:
            continue
        file = bz2.open(likelihood_timecourse_path,"rt")
        file.readline() # depth line!

        # total number of trajectories processed
        num_total = 0
        # total passed
        num_passed = 0

        annotated_mutations = []
        header = None
        printed_header = False
        loaded_avg_depths = False

        times_to_ignore = pt.samples_to_remove(population)


        for line in file:

            num_total += 1

            # parse timecourse data
            items = line.split(",")
            chromosome = items[0].strip()
            location = int(items[1])
            allele = items[2].strip()
            total_times = np.array([float(subitem) for subitem in items[3].split()])
            total_alts = np.array([float(subitem) for subitem in items[4].split()])

            total_depths = np.array([float(subitem) for subitem in items[5].split()])
            test_statistic = float(items[pvalue_idx].split()[0])
            pvalue = float(items[pvalue_idx].split()[1])

            # go back and figure out deletion and duplicatio
            # deletion_idx, fold_reduction, deletion_pvalue = tuple([float(subitem) for subitem in items[pvalue_idx+1].split()])
            #deletion_idx = long(deletion_idx)
            deletion_idx = len(total_times)
            deletion_pvalue=1
            fold_reduction=1

            #duplication_idx, fold_increase, duplication_pvalue = tuple([float(subitem) for subitem in items[pvalue_idx+2].split()])
            #duplication_idx = long(duplication_idx)
            duplication_idx=duplication_pvalue=fold_increase=1

            times = total_times[total_times<1000000]
            alts = total_alts[total_times<1000000]
            depths = total_depths[total_times<1000000]

            if times_to_ignore != None:
                times_to_ignore_idx=[np.where(times == t)[0] for t in times_to_ignore][0]

                times = np.delete(times, times_to_ignore_idx)
                alts = np.delete(alts, times_to_ignore_idx)
                depths = np.delete(depths, times_to_ignore_idx)

                total_times = np.delete(total_times, times_to_ignore_idx)
                total_alts = np.delete(total_alts, times_to_ignore_idx)
                total_depths = np.delete(total_depths, times_to_ignore_idx)


            # remove bad samples

            # load average depths if not leaded yet
            if not loaded_avg_depths:
                pop_avg_depths = depths
                loaded_avg_depths = True

            # create header if not created yet
            if header == None:
                print_strings = ['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 'Fold increase', 'Duplication pvalue', 'Passed?']
                for t in zip(total_times):
                    print_strings.append('AC:%d' % t)
                    print_strings.append('DP:%d' % t)

                header = ", ".join(print_strings)

            # annotate mutation
            gene_name, var_type = parse_file.annotate_variant(location, allele, gene_data, position_gene_map)

                # if pvalue is lower than threshold and not a weird insertion gene
            passed_str = 'FAIL'
            if (pvalue <= FDR) and (alts[0]*1.0/(depths[0]+(depths[0]==0)) < 0.2):
                #print(gene_name, var_type)

                # first estimate frequencies at good timepoints
                good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, deletion_idx, fold_reduction, deletion_pvalue)
                freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]

                masked_depth_ratios = depths[good_idxs]/pop_avg_depths[good_idxs]

                #interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)

                passed_str = 'PASS'
                num_passed += 1

            # print to CSV file
            print_strings = [str(location), gene_name, allele, var_type, str(test_statistic), str(pvalue), str(deletion_idx), str(fold_reduction), str(deletion_pvalue), str(duplication_idx), str(fold_increase), str(duplication_pvalue), passed_str]
            for alt,depth in zip(total_alts, total_depths):
                print_strings.append(str(alt))
                print_strings.append(str(depth))

            annotated_mutations.append((location, ", ".join(print_strings)))

        file.close()


        output_filename_template = pt.get_path() + "/data/timecourse_final/%s_annotated_timecourse.txt"
        # now write them out in sorted order
        output_file = open(output_filename_template % population,"w")
        output_file.write(header)
        output_file.write("\n")
        for location, mutation_str in sorted(annotated_mutations,key=lambda x: x[0]):
            output_file.write(mutation_str)
            output_file.write("\n")
        output_file.close()

        sys.stdout.write("%s: %d passed / %d\n" % (population, num_passed, num_total))

        total_num_passed += num_passed


# done iterating through populations
sys.stdout.write("%d total mutations!\n" % total_num_passed)



#sys.stdout.write("Filtering trajectories with FDR=%g\n" % FDR)
#input_filename_template = pt.get_path() + input_directory_prefix+"%s_likelihood_timecourse.txt" # change this later
#output_filename_template = pt.get_path() + output_directory_prefix+"%s_annotated_timecourse.txt"
