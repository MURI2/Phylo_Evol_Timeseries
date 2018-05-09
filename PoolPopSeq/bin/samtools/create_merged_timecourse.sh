#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -m n
#PBS -j oe

module load python/2.7.13

declare -a taxa=("B")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
#declare -a treats=("0")
#declare -a reps=("1")

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse_merged_fa

P_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.fna
D_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.fna
B_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna
C_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.fna
F_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna
J_fna=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.fna


for taxon in "${taxa[@]}"
do
  if [[ $taxon == "P" ]]; then
    REF_fna=$P_fna
    #continue
  elif [[ $taxon == "D" ]]
  then
    REF_fna=$D_fna
    #continue
  elif [[ $taxon == "B" ]]
  then
    REF_fna=$B_fna
    #continue
  elif [[ $taxon == "S" ]]
  then
    REF_fna=$B_fna
    #continue
  elif [[ $taxon == "C" ]]
  then
    REF_fna=$C_fna
    #continue
  elif [[ $taxon == "F" ]]
  then
    REF_fna=$F_fna
    #continue
  elif [[ $taxon == "J" ]]
  then
    REF_fna=$J_fna
    #continue
  else
    continue
  fi

  for treat in "${treats[@]}"
  do
    for rep in "${reps[@]}"
    do
      pop="Sample_L${treat}${taxon}${rep}"
      out_timecourse="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse_fa/${pop}_timecourse.txt"
      out_timecourse_merged="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse_merged_fa/${pop}_merged_timecourse.bz2"
      gd_files="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_fa/D"*"/${pop}-"*"/output/evidence/evidence.gd"
      cat ${out_timecourse} | python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/samtools/create_breseq_timecourse_no_repeat.py ${REF_fna} ${pop} ${gd_files} | bzip2 -c > ${out_timecourse_merged}
    done
  done
done
