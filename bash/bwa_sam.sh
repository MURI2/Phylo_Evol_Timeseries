#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# iRep and bPTR uses paired-end information, which breseq doesn't use
# so we'll need to map the reads using bwa and then clean them a little
# with samtools

module load bwa
module load samtools
module unload python
module load python/3.6.1

ref=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna

bwa index $ref
samtools faidx $ref

mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bwa_bam
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bwa_sam_merged
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/iRep
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bPTR
#declare -a strains=("B")
#declare -a treats=("0")
#declare -a reps=("2")
#declare -a times=("100")

declare -a strains=("S")
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
  if ls "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reads_clean_cutadapt/"*"${sample}"*"_R1_"*"_clean.fastq.gz" 1> /dev/null 2>&1; then
    for R1 in "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reads_clean_cutadapt/"*"${sample}"*"_R1_"*"_clean.fastq.gz"
    do
      R2="${R1/_R1_/_R2_}"
      bam_name="${R1/_R1/}"
      bam_name="${bam_name/fastq.gz/bam}"
      bam_name="${bam_name/reads_clean_cutadapt/bwa_bam}"
      bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
          | samtools sort -o $bam_name
    done

    bwa_sam_merged_out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bwa_sam_merged/Sample_${sample}.sam"
    samtools merge - "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bwa_bam/"*"_${sample}_"*".bam" \
        | samtools sort - | samtools view -h -o $bwa_sam_merged_out

    #/N/u/wrshoema/Carbonate/.local/bin/iRep -f $ref -s $bwa_sam_merged_out \
    #    -o "Sample_${sample}.iRep"

    /N/u/wrshoema/Carbonate/.local/bin/bPTR -f $ref -s $bwa_sam_merged_out \
        -m gc_skew -p 10000 -ff \
        -plot "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bPTR/Sample_${sample}.bPTR.pdf" \
        -o "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/bPTR/Sample_${sample}.bPTR.tsv"

  else
    continue

  fi
done
