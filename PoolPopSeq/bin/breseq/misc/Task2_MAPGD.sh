#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=8:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module rm gcc
module load gcc/4.9.2
module load gsl/1.15
module load samtools/0.1.19

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100

declare -a strains=("B" "C" "D" "F" "P")

for strain in "${strains[@]}"
do
  mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}"
  #echo "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}_names.txt"
  samtools merge "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_clean_paired_mapped_sort_NOdup_sort_${strain}_merged_merged.bam" \
    /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/Sample_L?$strain?/*clean_paired_mapped_sort_NOdup_sort_merged.bam
  samtools view -H "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_clean_paired_mapped_sort_NOdup_sort_${strain}_merged_merged.bam" \
    > "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_clean_paired_mapped_sort_NOdup_sort_${strain}_merged_merged.header"
  samtools mpileup -q 25 -Q 25 -B \
    /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/Sample_L?$strain?/*clean_paired_mapped_sort_NOdup_sort_merged.bam \
    | /N/dc2/projects/muri2/Task2/LTDE/MAPGD-master/bin/mapgd proview -H \
    "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_clean_paired_mapped_sort_NOdup_sort_${strain}_merged_merged.header" \
    | /N/dc2/projects/muri2/Task2/LTDE/MAPGD-master/bin/mapgd pool -a 22 -m 0.01 -o \
    "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_output_${strain}"

  OUT_Names="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/MAPGD_output/D100/MAPGD_output_${strain}/MAPGD_output_${strain}_merged_names.txt"
  if [ -s $OUT_Names ]
  then
    continue
  else
    echo /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/Sample_L?$strain?/*clean_paired_mapped_sort_NOdup_sort_merged.bam >> \
    $OUT_Names
  fi
done
