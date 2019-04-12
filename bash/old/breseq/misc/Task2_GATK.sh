#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load java


P=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.fna
D=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/Deinococcus_radiodurans_BAA816_genome.fa
B=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_168/AL009126.3.fa
C=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Caulobacter_crescentus_NA1000/CP001340.1.fa
F=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna


mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/GATK_output
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/GATK_output/D100

for bam in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/map_reads_trimmomatic/D100/*/*_clean_paired_mapped_sort_NOdup_sort_merged.bam
do
  line_folder="$(  echo "$bam" | cut -d"/" -f11-11)"
  OUTfolder="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/GATK_output/D100/${line_folder}"
  mkdir -p $OUTfolder

  bam_file="$(  echo "$bam" | cut -d"/" -f12-12)"
  NoExt="$(echo "${bam_file%.*}")"
  OUT="${OUTfolder}/${NoExt}"
  taxon="$(echo "$line_folder" | grep -Po ".(?=.{1}$)")"

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
  RefNoExt="$(echo "${REF%.*}")"
  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
      /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
      I=$bam \
      O="${OUT}_fixed.bam" \
      SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
      RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
      R=$REF \
      O="${RefNoExt}.dict"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
      /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
      -T RealignerTargetCreator \
      -R $REF \
      -I "${OUT}_fixed.bam" \
      -o "${OUT}_fixed.intervals"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
      -T IndelRealigner \
      -R $REF \
      -I "${OUT}_fixed.bam" \
      -targetIntervals "${OUT}_fixed.intervals" \
      --filter_bases_not_stored \
      -o "${OUT}_fixed_realigned.bam"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
      -T UnifiedGenotyper \
      -R $REF \
      -I "${OUT}_fixed_realigned.bam" \
      -glm BOTH -rf BadCigar \
      -o "${OUT}_fixed_realigned.vcf"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
      -T BaseRecalibrator -I "${OUT}_fixed_realigned.bam" \
      -R $REF \
      -rf BadCigar --filter_bases_not_stored -knownSites "${OUT}_fixed_realigned.vcf" \
      -o "${OUT}_fixed_realigned.recal_data.grp"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
      -T PrintReads -rf BadCigar \
      -R $REF \
      -I "${OUT}_fixed_realigned.bam" \
      -o "${OUT}_fixed_realigned_mapped.bam" \
      -BQSR "${OUT}_fixed_realigned.recal_data.grp"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
      -T UnifiedGenotyper \
      -R $REF \
      -I "${OUT}_fixed_realigned_mapped.bam" \
      -rf BadCigar \
      -o "${OUT}_fixed_realigned_mapped.vcf"

  java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -R $REF \
    -T VariantsToTable \
    -V "${OUT}_fixed_realigned_mapped.vcf" \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC \
    -o "${OUT}_fixed_realigned_mapped.txt"
done
