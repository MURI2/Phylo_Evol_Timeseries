#!/bin/bash

bash_rebreseq_scripts="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/rebreseq_scripts"
data_rebreseq="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq"
data_rebreseq_out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq_out"
data_rebreseq_err="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq_err"

mkdir -p $bash_rebreseq_scripts
mkdir -p $data_rebreseq
mkdir -p $data_rebreseq_out
mkdir -p $data_rebreseq_err

ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff

declare -a strains=("B")
declare -a treats=("0")
declare -a reps=("2")


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
  reads="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/reads_clean_cutadapt/"*"_${sample}_"*"_clean.fastq.gz"
  if (( ${#reads[@]} )); then
    pop="$(echo "$sample" | cut -d "_" -f1-1)"
    pop_gd="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq_jc/merged/${pop}.gd"
    bash_out="${bash_rebreseq_scripts}/${sample}_rebreseq.sh"
    if [ -f $bash_out ]; then
      rm $bash_out
    fi
    OUT_rebreseq="${data_rebreseq}/${sample}"
    mkdir -p $OUT_rebreseq
    OUT_rebreseq_out="${data_rebreseq_out}/${sample}.out"
    if [ -f $OUT_rebreseq_out ]; then
      rm $OUT_rebreseq_out
    fi
    OUT_rebreseq_err="${data_rebreseq_err}/${sample}.err"
    if [ -f $OUT_rebreseq_err ]; then
      rm $OUT_rebreseq_err
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#PBS -k o' >> $bash_out
    echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00' >> $bash_out
    echo '#PBS -M wrshoema@iu.edu' >> $bash_out
    echo '#PBS -m abe' >> $bash_out
    echo '#PBS -j oe' >> $bash_out
    echo '' >> $bash_out
    echo 'module load breseq' >> $bash_out
    # > data/breseq_output/${sample_name}.out 2> data/breseq_output/${sample_name}.err
    echo "breseq -j 8 -p --user-evidence-gd ${pop_gd} -o ${OUT_rebreseq} -r ${ref} ${reads} > ${OUT_rebreseq_out} 2> ${OUT_rebreseq_err}" >> $bash_out
    echo '' >> $bash_out

    qsub $bash_out


  else
    continue
  fi

done
