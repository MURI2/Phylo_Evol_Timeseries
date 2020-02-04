#!/bin/bash


create_timecourse=/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/create_timecourse.py


mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_merged
mkdir -p /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/timecourse_scripts

declare -a strains=("B")
#declare -a treats=("0" "1" "2")
#declare -a reps=("1" "2" "3" "4" "5")

declare -a treats=("0")
declare -a reps=("1")

declare -a pops=()

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      pops+=("${treat}${strain}${rep}")
    done
  done
done


for pop in "${pops[@]}"
do
  bash_out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/timecourse_scripts/${pop}_timecourse.sh"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi

  declare -a times=()
  for pop_sample in "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${pop}_"*"/data/reference.bam"
  do
    time="$(echo "$pop_sample" | cut -d "/" -f10-10 | cut -d "_" -f2-2)"
    times+=("${time}")
  done

  bam_files="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${pop}_"*"/data/reference.bam"
  ref="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/rebreseq/${pop}_100/data/reference.fasta"
  out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_merged/${pop}.pileup"
  if [ -f $out ]; then
    rm $out
  fi
  out_timecourse="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_merged/${pop}_timecourse.txt"

  echo '#!/bin/bash' >> $bash_out
  echo '#PBS -k o' >> $bash_out
  echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=12:00:00' >> $bash_out
  echo '#PBS -M wrshoema@iu.edu' >> $bash_out
  echo '#PBS -m abe' >> $bash_out
  echo '#PBS -j oe' >> $bash_out
  echo '' >> $bash_out
  echo 'module load samtools' >> $bash_out
  echo 'module load python' >> $bash_out
  echo "samtools mpileup -q10 -f ${ref} ${bam_files} > ${out}" >> $bash_out
  echo "cat ${out} | python ${create_timecourse} ${pop} ${times[@]} > $out_timecourse" >> $bash_out

  echo "${times[@]}"
  #qsub $bash_out
done
