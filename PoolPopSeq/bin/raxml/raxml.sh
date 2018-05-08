#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=8:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load raxml/8.0.26

cd /N/dc2/projects/muri2/Task2/PoolPopSeq/data/tree/

#raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# autoMRE \
#        -o Methanosarcina \
#        -s /N/dc2/projects/muri2/Task2/PoolPopSeq/data/tree/Task2_16S_rRNA.afa -n T20 \
#        -w /N/dc2/projects/muri2/Task2/PoolPopSeq/data/tree/

raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# autoMRE \
        -o Methanosarcina -n T20 \
        -s /Users/WRShoemaker/GitHub/Task2/PoolPopSeq/data/tree/Task2_16S_rRNA.afa \
        -w /Users/WRShoemaker/GitHub/Task2/PoolPopSeq/data/tree
