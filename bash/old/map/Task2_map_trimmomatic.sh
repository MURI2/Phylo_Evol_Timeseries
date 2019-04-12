#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load bwa/0.7.2
module load samtools/0.1.19

#J=
P=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.fna
D=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/Deinococcus_radiodurans_BAA816_genome.fa
B=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_168/AL009126.3.fa
C=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Caulobacter_crescentus_NA1000/CP001340.1.fa
F=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna

#bwa index $J
#samtools faidx $J
bwa index $P
samtools faidx $P
bwa index $D
samtools faidx $D
bwa index $B
samtools faidx $B
bwa index $C
samtools faidx $C
bwa index $F
samtools faidx $F

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100"

for folder in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/*/
do
  declare -a ARRAYreps=()
  declare -a ARRAYlines=()

  for file in $folder*_paired.fastq.gz
  do
    rep="$(  echo "$file" | cut -d'.' -f 1 | rev | cut -d"_" -f3 | rev)"
    ARRAYreps=("${ARRAYreps[@]}" "$rep")
    line="$(  echo "$file" | cut -d"_" -f2-6 | cut -d"/" -f3-5)"
    ARRAYlines=("${ARRAYlines[@]}" "$line")
  done
  REMDUPrep=($(printf "%s\n" "${ARRAYreps[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))
  REMDUPline=($(printf "%s\n" "${ARRAYlines[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))
  line_folder="$(  echo "$REMDUPline" | cut -d"/" -f1-1)"
  mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${line_folder}"
  taxon="$(echo "$folder" | grep -Po ".(?=.{2}$)")"

  if [[ $taxon == "P" ]]; then
    REF=$P
  elif [[ $taxon == "D" ]]
  then
    REF=$D
  elif [[ $taxon == "B" ]]
  then
    REF=$B
  elif [[ $taxon == "C" ]]
  then
    REF=$C
  elif [[ $taxon == "F" ]]
  then
    REF=$F
  else
    continue
  fi
  for rep in "${REMDUPrep[@]}"
  do
    InR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/${REMDUPline}_R1_${rep}_clean_paired.fastq.gz"
    InR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/${REMDUPline}_R2_${rep}_clean_paired.fastq.gz"
    Out="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${REMDUPline}_${rep}_clean_paired"
    bwa mem -t 4 $REF $InR1 $InR2 > "${Out}_mapped.sam"
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
  done
  OUTMerged="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${REMDUPline}_clean_paired_mapped_sort_NOdup_sort_merged"
  if [[ ${#REMDUPrep[@]} > 1 ]]; then
    samtools merge -h "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${REMDUPline}_001_clean_paired_mapped_sort_NOdup_sort.sam" \
      "${OUTMerged}.bam" \
      "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${REMDUPline}"_00?_clean_paired_mapped_sort_NOdup_sort.bam
  else
    cp "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/${REMDUPline}_001_clean_paired_mapped_sort_NOdup_sort.bam" \
      "${OUTMerged}.bam"
  fi
  samtools index "${OUTMerged}.bam"
  samtools view -h -o "${OUTMerged}.sam" "${OUTMerged}.bam"
done
