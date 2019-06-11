#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=4:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

# version 8.2.11
module load raxml

raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o NC_005042.1.353331-354795 -# autoMRE \
    -s /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/tree/arb-silva.de_2019-06-07_id664129.fasta \
    -n Task2 -w /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/tree



# -T = number of threads
# -f = specifies bootstrapping algorithm with ML generating tree at same time
# -m = substitution model, generalized time reversible gamma
# -p = starts tree randomly
# -x = starts tree randomly
# -o = outgroup (name after fasta entry)
# -#  = determines number of bootstrap replicates
# -s = aligned fasta file to be analyzed
# -n = name of output file
