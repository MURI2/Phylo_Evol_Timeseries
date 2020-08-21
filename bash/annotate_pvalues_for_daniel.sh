#!/bin/bash


declare -a strains=("B")
declare -a treats=("0")
declare -a reps=("1")

#declare -a treats=("0")
#declare -a reps=("1")

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
  bash_out="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/bash/annotate_pvalues_scripts/${pop}_annotate_pvalues.sh"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi

  depth_timecourse_filename="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_depth/${pop}_depth_timecourse.bz"
  snp_timecourse_filename="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_snp/${pop}_snp_timecourse.bz"
  indel_timecourse_filename="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_indel/${pop}_indel_timecourse.bz"
  likelihood_timecourse_filename="/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/timecourse_likelihood/${pop}_likelihood_timecourse.bz"

  echo '#!/bin/bash' >> $bash_out
  echo '#PBS -k o' >> $bash_out
  echo '#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=12:00:00' >> $bash_out
  echo '#PBS -M wrshoema@iu.edu' >> $bash_out
  echo '#PBS -m abe' >> $bash_out
  echo '#PBS -j oe' >> $bash_out
  echo '' >> $bash_out
  echo 'module unload python' >> $bash_out
  echo 'module load python/3.6.1' >> $bash_out

  echo "python /N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/Python/annotate_pvalues.py ${depth_timecourse_filename} ${snp_timecourse_filename} ${indel_timecourse_filename} ${likelihood_timecourse_filename}" >> $bash_out

  qsub $bash_out
done
