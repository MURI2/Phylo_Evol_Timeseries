#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load cutadapt

declare -a days=("100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")
#declare -a days=("1000")

#A, F, J, P

# trim data and remove adaptors.
mkdir -p "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reads_clean_cutadapt"

for day in "${days[@]}"
do
  for R1 in "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/illumina_runs_rename/"*"_"*"P"*"_${day}_"*"_R1_"*".fastq.gz";
  do
    if [ ! -f $R1 ]; then
      continue
    fi

    R2="${R1/_R1_/_R2_}"
    R1_clean="${R1/.fastq.gz/_clean.fastq.gz}"
    R1_clean="${R1_clean/illumina_runs_rename/reads_clean_cutadapt}"
    R2_clean="${R2/.fastq.gz/_clean.fastq.gz}"
    R2_clean="${R2_clean/illumina_runs_rename/reads_clean_cutadapt}"

    cutadapt -q 30,30 -u 20 \
              -a file:/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/transposase.fa \
              --minimum-length 20 \
              -o $R1_clean \
              -p $R2_clean \
              $R1 $R2
  done
done
