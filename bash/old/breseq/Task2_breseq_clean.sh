#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=48:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load gcc/4.9.4
#module load samtools
#module load python


DAY=D100
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}"

#for sample in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/D100/*;
#do
#  line="$(echo "$sample" | cut -d "/" -f11-11)"
#  OUT="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/D100/${line}"
#  #mkdir -p $OUT
#  #cp "${sample}/output/evidence/annotated.gd" "${OUT}/annotated.gd"
#  #cp "${sample}/output/evidence/evidence.gd" "${OUT}/evidence.gd"
#  #cp "${sample}/output/output.gd" "${OUT}/output.gd"
#  #samtools mpileup "${sample}/data/reference.bam" > "${sample}/data/coverage.txt"
#  #python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -c -i \
#  #  "${sample}/data/coverage.txt" -o "${sample}/data/coverage_merged.txt"
#  #"${OUT}/evidence.gd"
#done

#python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/calculateNum.py
#python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/test.py

for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f11-11)"
  bash_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/breseq/line_scripts_gbk/${line}_breseq.sh"
  if [ ! -f $bash_out ]; then
    OUT="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}/${line}"
    OUT_essentials="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}/${line}"
    mkdir -p $OUT_essentials

    OUT_annotated="${OUT}/output/evidence/annotated.gd"
    essentials_annotated="${OUT_essentials}/annotated.gd"
    cp $OUT_annotated $essentials_annotated

    OUT_evidence="${OUT}/output/evidence/evidence.gd"
    essentials_evidence="${OUT_essentials}/evidence.gd"
    cp $OUT_evidence $essentials_evidence

    OUT_output="${OUT}/output/output.gd"
    essentials_output="${OUT_essentials}/output.gd"
    cp $OUT_output $essentials_output

    bam="${OUT}/data/reference.bam"
    coverage="${OUT}/data/coverage.txt"
    if [ -fe $bam ]; then
        rm $bam
    fi
    samtools mpileup $bam > $coverage

    coverage_clean="${OUT_essentials}/coverage_clean.txt"
    python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -c -i $coverage -o $coverage_clean
  fi
done
