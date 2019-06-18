#!/bin/bash

mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq_json

declare -a strains=("B" "S")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")

declare -a times=("100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")

declare -a samples=()

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      for time in "${times[@]}"
      do
        samples+=("${treat}${strain}${rep}_${time}")
      done
    done
  done
done


for sample in "${samples[@]}"
do
  json="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${sample}/data/summary.json"
  new_json="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq_json/${sample}.json"
  if [ ! -f $json ]; then
    continue
  fi
  cp $json $new_json
done
