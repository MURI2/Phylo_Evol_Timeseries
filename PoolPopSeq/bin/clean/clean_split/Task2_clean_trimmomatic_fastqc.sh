#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=6:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc

Sample_date=D200

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic_quality/"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic_quality/${Sample_date}/"


for file in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/"*/*.fastq.gz
do
  sample="$(echo "$file" | cut -d "/" -f11-11)"
  OUTdir="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic_quality/${Sample_date}/${sample}"
  mkdir -p $OUTdir
  #echo "$file"
  fastqc "$file" --outdir=$OUTdir
done
