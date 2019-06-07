'''
This script was originally written by B.H. Good for the publication
The dynamics of molecular evolution over 60,000 generations
Nature volume 551, pages 45-50, doi:10.1038/nature24287
This script has been modified with permission and is free for use
under a GNU General Public License v2.0
'''


# This script takes the output from our second breseq run
# (the "rebreseq" step") and calculates a timecourse for
# each junction. Because breseq will call junctions slightly
# differently at different timepoints, this script includes
# a fuzzy matching algorithm for merging junction candidates
# that are likely to be the same.

# the idea is that allele reads get added together and references get averaged

import numpy
import sys
#import population_parameters
from operator import itemgetter, attrgetter, methodcaller
from parse_file import parse_repeat_list, get_repeat_idx, get_closest_repeat_idx
from math import fabs

# File containing the reference fasta file
reference_filename = sys.argv[1]
# File containing the reference gbk/gbff file
reference_gbk_filename = sys.argv[2]
# The population to compile timecourses for
population = sys.argv[3]
# The output (gd) files from breseq
gd_filenames = sys.argv[4:]


# How far can two simple indels be before they are merged
INDEL_EDGE_TOLERANCE = 0
# How long can a simple indel be before it is treated differently?
INDEL_LENGTH_TOLERANCE = 100
# How far can two similar IS junctions be before they are merged
REPEAT_EDGE_TOLERANCE = 20

# How far can two similar inversion elements be before they are merged
OTHER_EDGE_TOLERANCE = 20


INDEL = 0
REPEAT = 1
OTHER = 2

# Construct the list of samples for this population
#full_sample_list = population_parameters.sample_names[population]
#full_sample_times = population_parameters.sample_times[population]

sample_list = []
gd_filename_map = {}
#parse sample names from filenames
for gd_filename in gd_filenames:
    sample_name = gd_filename.split("/")[9].strip()
    if sample_name.startswith(population):
        sample_name = sample_name.split("_",1)[1]
    sample_list.append(sample_name)
    gd_filename_map[sample_name] = gd_filename

print(gd_filename_map)

sample_list = sorted(sample_list)
