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


import json



json_path = pt.get_path() + '/data/rebreseq_json/'


coverages_all = []
for filename in os.listdir(json_path):

    if filename.endswith(".json"):


        filepath = '%s%s' % (json_path, filename)

        with open(filepath) as f:
            data = json.load(f)

        #print(data.keys())

        #print()

        contigs = data['references']['reference'].keys()

        coverages = [data['references']['reference'][contig]['coverage_average'] for contig in contigs]


        coverages_all.extend(coverages)



        # coverage_average


mean_coverage = np.mean(coverages_all)

se_coverage = np.std(coverages_all) / np.sqrt(len(coverages_all))


print(mean_coverage, se_coverage)
