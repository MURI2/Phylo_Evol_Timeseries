declare -a strains=("A")
declare -a treats=("2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300" "400" "500" "600" "700")

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
  echo $sample
  samtools depth "/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries/data/breseq/${sample}/data/reference.bam"  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
done
