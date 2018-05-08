#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load python
module load gcc
module load cutadapt
module load java
module load fastqc

#-g CAAGCAGAAGACGGCATACGA
#-g AATGATACGGCGACCACCGA

Sample_date=D100

# raw data cleaned of adapters
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/${Sample_date}"

for folder in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/"*/
do
  declare -a ARRAYlines=()

  for file in $folder*.fastq.gz
  do
    if [[ $file == *"_R1_"* ]]; then
      line="$(  echo "$file" | cut -d"/" -f11-13)"
      ARRAYlines=("${ARRAYlines[@]}" "$line")
    fi
  done
  line_folder="$(  echo "$folder" | cut -d"/" -f11-12)"
  mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/${Sample_date}/${line_folder}"
  for ARRAYline in "${ARRAYlines[@]}"
  do
    R1=$ARRAYline
    R2="${R1/_R1_/_R2_}"
    InR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/${R1}"
    InR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/${R2}"
    OutR1Name="${R1/.fastq.gz/_clean.fastq.gz}"
    OutR2Name="${R2/.fastq.gz/_clean.fastq.gz}"
    OutR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/${Sample_date}/${OutR1Name}"
    OutR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/${Sample_date}/${OutR2Name}"
    adaptor1="$(  echo "$OutR1" | cut -d"_" -f4-4 | cut -d"-" -f1-1)"
    adaptor2="$(  echo "$OutR1" | cut -d"_" -f4-4 | cut -d"-" -f2-2)"
    #   -b AGATCGGAAGAGC -B AGATCGGAAGAGC -b $adaptor1 -B $adaptor2 \
    # -b AATGATACGGCGACCACCGA -B CAAGCAGAAGACGGCATACGA \
    #-b "CAAGCAGAAGACGGCATACGAGAT${adaptor1}GTCTCGTGGGCTCGG" \
    #-B "AATGATACGGCGACCACCGAGATCTACAC${adaptor2}TCGTCGGCAGCGTC" \
    #-b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -B GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    cutadapt -q 30,30 -u 10 \
      -b $adaptor1 -B $adaptor1 -b $adaptor2 -B $adaptor2 \
      -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -B GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
      -o $OutR1 -p $OutR2 $InR1 $InR2
  done
done

# fastqc of cleaned data
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_quality/
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_quality/${Sample_date}/"

for file in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/${Sample_date}/"*/*.fastq.gz
do
  sample="$(echo "$file" | cut -d "/" -f11-11)"
  OUTdir="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_quality/${Sample_date}/${sample}"
  mkdir -p $OUTdir
  fastqc "$file" --outdir=$OUTdir
done
