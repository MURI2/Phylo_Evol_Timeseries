#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=5:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load gcc/4.9.4
module load samtools/1.3.1
module load python/2.7.3

DAY=D200

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage/${DAY}"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage_clean/${DAY}"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage_clean_figs/${DAY}"

for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${DAY}/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f11-11)"
  taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
  OUT="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}/${line}"
  bam="${OUT}/data/reference.bam"
  OUT_cov="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage/${DAY}/${line}_coverage.txt"
  samtools mpileup ${bam} > ${OUT_cov}
  # next, clean coverage and take average and std dev every 100 bases
  OUT_cov_clean="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage_clean/${DAY}/${line}_coverage_clean.txt"
  OUT_cov_clean_plt="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/coverage_clean_figs/${DAY}/${line}_coverage_clean"
  python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -m -i ${OUT_cov} -o ${OUT_cov_clean}
done
