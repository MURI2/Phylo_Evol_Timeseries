#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=5:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/illumina_runs_rename/

while read -r first second; do
    old="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/illumina_runs/${first}"
    new="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/illumina_runs_rename/${second}"
    #echo $old
    #echo $new
    cp $old $new
done < /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/new_sample_names.txt
