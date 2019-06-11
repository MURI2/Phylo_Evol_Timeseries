#!/bin/bash

module load breseq

#declare -a strains=("B" "C" "D" "F" "J" "P" "S")
#declare -a treats=("0" "1" "2")
#declare -a reps=("1" "2" "3" "4" "5")

declare -a strains=("F" "J" "P")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")



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


mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/all
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/merged


for pop in "${pops[@]}"
do
  for time in "${times[@]}"
  do
    evidence_time="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq/${pop}_${time}/output/evidence/evidence.gd"
    if [ -f $evidence_time ] ; then
      junction_output="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/all/${pop}_${time}.gd"
      cat $evidence_time | grep 'JC\|#' > $junction_output
    fi
  done
  merged_junction_output="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/merged/${pop}.gd"
  if [ -f $junction_merged_output ] ; then
    rm $junction_merged_output
  fi
  evidence="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/all/${pop}"*".gd"
  gdtools UNION -o $merged_junction_output -e $evidence
done
