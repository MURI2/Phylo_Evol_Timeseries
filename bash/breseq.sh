#!/bin/bash

bin_breseq_scripts="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/breseq_scripts"
data_breseq="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq"
data_breseq_out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_out"
data_breseq_err="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_err"

mkdir -p $bin_breseq_scripts
mkdir -p $data_breseq
mkdir -p $data_breseq_out
mkdir -p $data_breseq_err

A_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Arthrobacter_sp_KBS0703/GCF_002008315.2_ASM200831v2_genomic.gbff
B_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.gbff
C_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCF_000022005.1_ASM2200v1_genomic.gbff
D_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCF_000008565.1_ASM856v1_genomic.gbff
F_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Pedobacter_sp_KBS0701/GCF_005938645.2_ASM593864v2_genomic.gbff
J_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/GCF_005937955.2_ASM593795v2_genomic.gbff
P_gbk=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/GCF_005938045.2_ASM593804v2_genomic.gbff

#declare -a strains=("A" "B" "C" "D" "F" "J" "P" "S")
# Done: B,S,C,D,
declare -a treats=("0")
declare -a strains=("F")

declare -a reps=("3")

#declare -a times=("3000")
declare -a times=("100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")


#declare -a samples=("1C4_100" "1D3_100" "1F2_100" "1F3_100" "1F5_100" "1P3_100" "2B3_100" "2C4_100" "2C5_100" "2D4_100" "2F5_100" "2P5_100")
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
  reads="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reads_clean_cutadapt/"*"_${sample}_"*"_clean.fastq.gz"
  files_test=( $reads )
  #if (( ${#reads[@]} )); then
  if [ -e "${files_test[0]}" ]; then
    bash_out="${bin_breseq_scripts}/${sample}_breseq.sh"
    if [ -f $bash_out ]; then
      rm $bash_out
    fi
    OUT_breseq="${data_breseq}/${sample}"
    # delete breseq directory if it exists
    if [ -d $OUT_breseq ]
    then
      rm -r $OUT_breseq
    fi
    mkdir -p $OUT_breseq

    OUT_breseq_out="${data_breseq_out}/${sample}.out"
    if [ -f $OUT_breseq_out ]; then
      rm $OUT_breseq_out
    fi
    OUT_breseq_err="${data_breseq_err}/${sample}.err"
    if [ -f $OUT_breseq_err ]; then
      rm $OUT_breseq_err
    fi

    # get reference
    if [[ $sample == *"B"* ]]; then
      gbk=$B_gbk
    elif [[ $sample == *"S"* ]]; then
      gbk=$B_gbk
    elif [[ $sample == *"C"* ]]; then
      gbk=$C_gbk
    elif [[ $sample == *"D"* ]]; then
      gbk=$D_gbk
    elif [[ $sample == *"F"* ]]; then
      gbk=$F_gbk
    elif [[ $sample == *"J"* ]]; then
      gbk=$J_gbk
    elif [[ $sample == *"P"* ]]; then
      gbk=$P_gbk
    elif [[ $sample == *"A"* ]]; then
      gbk=$A_gbk
    else
      continue
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#PBS -k o' >> $bash_out
    echo '#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=10:00:00' >> $bash_out
    #echo '#PBS -M wrshoema@iu.edu' >> $bash_out
    echo '#PBS -m abe' >> $bash_out
    echo '#PBS -j oe' >> $bash_out
    echo '' >> $bash_out
    echo 'module load breseq' >> $bash_out
    echo "breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o ${OUT_breseq} -r ${gbk} ${reads} > ${OUT_breseq_out} 2> ${OUT_breseq_err}" >> $bash_out
    echo '' >> $bash_out

    qsub $bash_out

  else
    continue
  fi

done
