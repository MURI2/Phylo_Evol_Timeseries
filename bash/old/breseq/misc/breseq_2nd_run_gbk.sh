#!/bin/bash

#DAY=D100
DAY=$1

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_2nd_run_gbk
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_2nd_run_gbk/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk_out
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk_out/${DAY}"

P_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.gbk
D_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.gbff
B_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff
C_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.gbff
F_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.gbk
J_gbk=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.gbk


for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${DAY}/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f11-11)"
  pop="$(echo "$line_path" | cut -d "/" -f11-11 | cut -d "-" -f1-1)"
  pop_gd="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_1st_run_gbk_jc/merged/${pop}.gd"
  bash_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_2nd_run_gbk/${DAY}/${line}_breseq.sh"
  #if [ -f $bash_out ]; then
  #  rm $bash_out
  #fi
  taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
  if [[ $taxon == "P" ]]; then
    REF_gbk=$P_gbk
    #continue
  elif [[ $taxon == "D" ]]
  then
    REF_gbk=$D_gbk
    #continue
  elif [[ $taxon == "B" ]]
  then
    REF_gbk=$B_gbk
    #continue
  elif [[ $taxon == "S" ]]
  then
    #REF_gbk=$B_gbk
    continue
  elif [[ $taxon == "C" ]]
  then
    REF_gbk=$C_gbk
    #continue
  elif [[ $taxon == "F" ]]
  then
    REF_gbk=$F_gbk
    #continue
  elif [[ $taxon == "J" ]]
  then
    REF_gbk=$J_gbk
    #continue
  else
    continue
  fi

  reads="${line_path}/"*_clean_paired.fastq.gz
  OUT_gbk="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk/${DAY}/${line}"
  mkdir -p $OUT_gbk
  OUT_folder="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk_out/${DAY}/${line}"
  mkdir -p $OUT_folder
  OUT_stdout="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk_out/${DAY}/${line}/${line}.out"
  OUT_stderr="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_2nd_run_gbk_out/${DAY}/${line}/${line}.err"

  echo '#!/bin/bash' >> $bash_out
  echo '#PBS -k o' >> $bash_out
  echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00' >> $bash_out
  echo '#PBS -M wrshoema@umail.iu.edu' >> $bash_out
  #echo '#PBS -m abe' >> $bash_out
  echo '#PBS -m n' >> $bash_out
  echo '#PBS -j oe' >> $bash_out
  echo '' >> $bash_out
  #echo 'module load python/2.7.3' >> $bash_out
  #echo 'module unload gcc/4.9.4' >> $bash_out
  #echo 'module load gcc/4.9.2' >> $bash_out
  #echo 'module load bowtie2/2.2.6' >> $bash_out
  #echo 'module load r/3.3.1' >> $bash_out
  echo 'module load breseq' >> $bash_out
  echo "breseq -j 8 -p --user-evidence-gd ${pop_gd} -o ${OUT_gbk} -r ${REF_gbk} ${reads} > ${OUT_stdout} 2> ${OUT_stderr}" >> $bash_out
  echo '' >> $bash_out

  qsub $bash_out
done
