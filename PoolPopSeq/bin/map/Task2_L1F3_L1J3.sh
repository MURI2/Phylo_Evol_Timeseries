#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load python
module load gcc
module load cutadapt
module load bwa/0.7.2
module load samtools/0.1.19
module load vcftools

#-g CAAGCAGAAGACGGCATACGA
#-g AATGATACGGCGACCACCGA

OutR1Paired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/Sample_L1J3/L1J3_CGGAGCCT-GTAAGGAG_L002_R1_001_clean_paired.fastq.gz"
OutR2Paired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/Sample_L1J3/L1J3_CGGAGCCT-GTAAGGAG_L002_R2_001_clean_paired.fastq.gz"


REF=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna
bwa index $REF
samtools faidx $REF

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/Sample_L1J3"
Out="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/Sample_L1J3/L1J3_CGGAGCCT-GTAAGGAG_L002_001_clean_paired"

bwa mem -t 4 $REF $OutR1Paired $OutR2Paired > "${Out}_mapped.sam"
# mapped reads
samtools view -F 4 -bT $REF "${Out}_mapped.sam" > "${Out}_mapped.bam"
# unmapped reads
samtools view -f 4 -bT $REF "${Out}_mapped.sam" > "${Out}_unmapped.bam"

samtools sort "${Out}_mapped.bam" "${Out}_mapped_sort"
samtools index "${Out}_mapped_sort.bam"
samtools rmdup "${Out}_mapped_sort.bam" "${Out}_mapped_sort_NOdup.bam"
samtools index "${Out}_mapped_sort_NOdup.bam"
samtools sort "${Out}_mapped_sort_NOdup.bam" "${Out}_mapped_sort_NOdup_sort"
samtools index "${Out}_mapped_sort_NOdup_sort.bam"
samtools view -h -o "${Out}_mapped_sort_NOdup_sort.sam" "${Out}_mapped_sort_NOdup_sort.bam"
# same thing for unmapped reads
samtools sort "${Out}_unmapped.bam" "${Out}_unmapped_sort"
samtools index "${Out}_unmapped_sort.bam"
samtools rmdup "${Out}_unmapped_sort.bam" "${Out}_unmapped_sort_NOdup.bam"
samtools index "${Out}_unmapped_sort_NOdup.bam"
samtools sort "${Out}_unmapped_sort_NOdup.bam" "${Out}_unmapped_sort_NOdup_sort"
samtools index "${Out}_unmapped_sort_NOdup_sort.bam"
