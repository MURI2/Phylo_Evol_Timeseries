#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=6:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic_fastqc
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic_fastqc/D100/

for file in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/*/*_clean_paired_mapped_sort_NOdup_sort.bam
do
  sample="$(echo "$file" | cut -d "/" -f11-11)"
  OUTdir="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic_fastqc/D100/${sample}"
  mkdir -p $OUTdir
  #echo $sample
  fastqc "$file" --outdir=$OUTdir
done
