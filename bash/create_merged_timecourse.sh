#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=15:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -m n
#PBS -j oe

module load python

declare -a strains=("F" "J" "P" "S")
declare -a treats=("0" "1")
declare -a reps=("1" "2" "3" "4" "5")

declare -a pops=()

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      pops+=("${treat}${strain}${rep}")
    done
  done
done

create_breseq_timecourse=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/create_breseq_timecourse.py

for pop in "${pops[@]}"
do
  gd_files="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${pop}_"*"/output/evidence/evidence.gd"
  #ref="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${pop}_100/data/reference.fasta"
  cat "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_merged/${pop}_timecourse.txt" | python $create_breseq_timecourse $pop $gd_files | bzip2 -c > "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_merged/${pop}_merged_timecourse.bz"
done
