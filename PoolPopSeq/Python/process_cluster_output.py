#############################
#
# This script processes output from cluster to obtain final set
# of mutation trajectories that are used for downstream analysis
#
#############################

import os, imp, sys
from datetime import date

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")
parse_file = imp.load_source('parse_file', mydir + 'bin/samtools/parse_file.py')

storage_directory='cluster_output_files/'
final_directory = mydir + 'data/final_directory'

#treatments = ['0', '1', '2']
#reps = ['1', '2', '3', '4', '5']
treatments = ['0']
reps = ['1']
#strains = ['B', 'C', 'D', 'F', 'J', 'P', 'S']
strains = ['B']

def process_output():
    for strain in strains:
        for treatment in treatments:
            for rep in reps:
                sys.stdout.write("\nProcessing Sample_L%s%s%s...\n" % (treatment, strain, rep))
                merged_timecourse_filename = 'Sample_L%s%s%s_merged_timecourse.bz2' % (treatment, strain, rep)
                depth_timecourse_filename = 'Sample_L%s%s%s_depth_timecourse.bz2' % (treatment, strain, rep)
                snp_timecourse_filename = 'Sample_L%s%s%s_snp_timecourse.bz2' % (treatment, strain, rep)
                indel_timecourse_filename = 'Sample_L%s%s%s_indel_timecourse.bz2' % (treatment, strain, rep)
                likelihood_timecourse_filename = 'Sample_L%s%s%s_likelihood_timecourse.txt' % (treatment, strain, rep)

                merged_timecourse_path = mydir + 'data/timecourse_merged_fa/' + merged_timecourse_filename
                depth_timecourse_path = mydir + 'data/timecourse_depth/' + depth_timecourse_filename
                snp_timecourse_path = mydir + 'data/timecourse_snp/' + snp_timecourse_filename
                indel_timecourse_path = mydir + 'data/timecourse_indel/' + indel_timecourse_filename
                likelihood_timecourse_path = mydir + 'data/timecourse_likelihood/' + likelihood_timecourse_filename

                # Filter SNPs and calculate avg depth per sample
                sys.stdout.write('Filtering SNPS and calculating depth...\n')
                filter_snps_path = mydir + 'Python/filter_snps_and_calculate_depth.py'
                return_val = os.system('python %s %s %s %s' % (filter_snps_path, merged_timecourse_path, depth_timecourse_path, snp_timecourse_path))
                if return_val==0:
                    sys.stdout.write('Done!\n')
                else:
                    sys.stdout.write("Error!\n")

                # Call indels
                sys.stdout.write("Calling indels...\n")
                call_indels_path = mydir + 'Python/call_indels.py'
                return_val = os.system('python %s %s %s' % (call_indels_path, merged_timecourse_path, indel_timecourse_path))
                if return_val==0:
                    sys.stdout.write('Done!\n')
                else:
                    sys.stdout.write("Error!\n")

                # Annotate pvalues
                #sys.stdout.write("Calculating pvalues...\n")
                pvals_path = mydir + 'cpp/annotate_pvalues'
                print 'bzcat %s %s %s | %s > %s' % \
                    (depth_timecourse_path, snp_timecourse_path, indel_timecourse_path, \
                    pvals_path, likelihood_timecourse_path)
                return_val = os.system('bzcat %s %s %s | %s > %s' % \
                    (depth_timecourse_path, snp_timecourse_path, indel_timecourse_path, \
                    pvals_path, likelihood_timecourse_path))
                if return_val==0:
                    sys.stdout.write('Done!\n')
                else:
                    sys.stdout.write("Error!\n")

                #sys.stdout.write("\n\nTrajectory post-processing output for LTEE metagenomic sequencing project\n")
                #sys.stdout.write("Date: %s\n\n" % str(date.today()))
                #os.system('mkdir -p %s' % final_directory)

                # Filter and annotate timecourse
                #sys.stdout.write('Performing final filtering and annotation step...\n')
                #return_val = os.system('python combined_annotate_timecourse.py %s %s' % (storage_directory, final_directory))
                #if return_val==0:
                #    sys.stdout.write('Done!\n')
                #else:
                #    sys.stdout.write("Error!\n")

                # Infer trajectory states in well-mixed HMM
                #sys.stdout.write('Inferring trajectory states in well-mixed HMM...\n')
                #return_val = os.system('python calculate_well_mixed_hmm_wrapper.py')
                #if return_val==0:
                #    sys.stdout.write('Done!\n')
                #else:
                #    sys.stdout.write("Error!\n")



process_output()



# Infer trajectory states in clade HMM
#sys.stdout.write('Inferring trajectory states in clade HMM...\n')
#return_val = os.system('python calculate_clade_hmm_wrapper.py')
#if return_val==0:
#    sys.stdout.write('Done!\n')
#else:
#    sys.stdout.write("Error!\n")
