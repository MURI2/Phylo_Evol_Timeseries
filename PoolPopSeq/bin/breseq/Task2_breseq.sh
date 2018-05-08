#!/bin/bash

DAY=D100

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_scripts
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_scripts/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}"

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_fna
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_fna/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_fna
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_fna/${DAY}"


P_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.gbk
D_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.gbff
B_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff
C_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.gbff
F_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.gbk
J_gbk=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.gbk

P_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.fna
D_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.fna
B_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna
C_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.fna
F_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna
J_fna=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.fna


for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${DAY}/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f11-11)"
  bash_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq/breseq_scripts/${DAY}/${line}_breseq.sh"
  #echo $line
  if [ ! -f $bash_out ]; then
    taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
    if [[ $taxon == "P" ]]; then
      #REF_gbk=$P_gbk
      #REF_fna=$P_fna
      continue
    elif [[ $taxon == "D" ]]
    then
      #REF_gbk=$D_gbk
      #REF_fna=$D_fna
      continue
    elif [[ $taxon == "B" ]]
    then
      REF_gbk=$B_gbk
      REF_fna=$B_fna
      #continue
    elif [[ $taxon == "S" ]]
    then
      REF_gbk=$B_gbk
      REF_fna=$B_fna
      #continue
    elif [[ $taxon == "C" ]]
    then
      #REF_gbk=$C_gbk
      #REF_fna=$C_fna
      continue
    elif [[ $taxon == "F" ]]
    then
      #REF_gbk=$F_gbk
      #REF_fna=$F_fna
      continue
    elif [[ $taxon == "J" ]]
    then
      #REF_gbk=$J_gbk
      #REF_fna=$J_fna
      continue
    else
      continue
    fi

    reads="${line_path}/"*_clean_paired.fastq.gz
    OUT_gbk="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}/${line}"
    OUT_bwa_fna="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_fna/${DAY}/${line}"
    OUT_gbk_essentials="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}/${line}"
    mkdir -p $OUT_gbk
    mkdir -p $OUT_bwa_fna
    mkdir -p $OUT_gbk_essentials

    echo '#!/bin/bash' >> $bash_out
    echo '#PBS -k o' >> $bash_out
    echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00' >> $bash_out
    echo '#PBS -M wrshoema@umail.iu.edu' >> $bash_out
    echo '#PBS -m abe' >> $bash_out
    echo '#PBS -j oe' >> $bash_out
    echo '' >> $bash_out
    echo 'module load python' >> $bash_out
    echo 'module load gcc/4.9.4' >> $bash_out
    echo 'module load bowtie2/2.2.6' >> $bash_out
    echo 'module load intel' >> $bash_out
    echo 'module load curl' >> $bash_out
    echo 'module load java' >> $bash_out
    echo 'module load R/3.3.1' >> $bash_out
    echo 'module load breseq/0.27' >> $bash_out
    echo 'module load bwa/0.7.2' >> $bash_out
    echo 'module load samtools/1.3.1' >> $bash_out
    echo '' >> $bash_out
    #echo "breseq -j 8 -p -o ${OUT_gbk} -r ${REF_gbk} ${reads}" >> $bash_out
    # run bwa on each set of files (R1 and R2)
    # iterate through each R1 and R2
    echo '' >> $bash_out
    echo "bwa index ${REF_fna}" >> $bash_out
    echo "samtools faidx ${REF_fna}" >> $bash_out
    for read in "${line_path}/"*_R1_*_clean_paired.fastq.gz
    do
      if [[ $read == *"_R1_"* ]]; then
        R1=$read
        R2="${read/_R1_/_R2_}"
        line_rep="$(  echo "$read" | cut -d"/" -f12-13 | cut -d"_" -f1-6)"
        line_rep="${line_rep/_R1_/_}"
        out_line_rep="${OUT_bwa_fna}/Sample_${line_rep}"
        echo "bwa mem -t 4 ${REF_fna} ${R1} ${R2} > "${out_line_rep}_mapped.sam"" >> $bash_out
        # mapped reads
        echo "samtools view -F 4 -bT ${REF_fna} "${out_line_rep}_mapped.sam" > "${out_line_rep}_mapped.bam"" >> $bash_out
        # unmapped reads
        echo "samtools view -f 4 -bT ${REF_fna} "${out_line_rep}_mapped.sam" > "${out_line_rep}_unmapped.bam"" >> $bash_out
        echo "samtools sort "${out_line_rep}_mapped.bam" -o "${out_line_rep}_mapped_sort.bam""  >> $bash_out
        echo "samtools index "${out_line_rep}_mapped_sort.bam"" >> $bash_out
        echo "samtools rmdup "${out_line_rep}_mapped_sort.bam" "${out_line_rep}_mapped_sort_NOdup.bam"" >> $bash_out
        echo "samtools index "${out_line_rep}_mapped_sort_NOdup.bam"" >> $bash_out
        echo "samtools sort "${out_line_rep}_mapped_sort_NOdup.bam" -o "${out_line_rep}_mapped_sort_NOdup_sort.bam"" >> $bash_out
        echo "samtools index "${out_line_rep}_mapped_sort_NOdup_sort.bam"" >> $bash_out
        # same thing for unmapped reads
        echo "samtools sort "${out_line_rep}_unmapped.bam" -o "${out_line_rep}_unmapped_sort.bam"" >> $bash_out
        echo "samtools index "${out_line_rep}_unmapped_sort.bam"" >> $bash_out
        echo "samtools rmdup "${out_line_rep}_unmapped_sort.bam" "${out_line_rep}_unmapped_sort_NOdup.bam"" >> $bash_out
        echo "samtools index "${out_line_rep}_unmapped_sort_NOdup.bam"" >> $bash_out
        echo "samtools sort "${out_line_rep}_unmapped_sort_NOdup.bam" -o "${out_line_rep}_unmapped_sort_NOdup_sort.bam"" >> $bash_out
        echo "samtools index "${out_line_rep}_unmapped_sort_NOdup_sort.bam"" >> $bash_out
      fi
    done
    echo '' >> $bash_out
    # merge the bams
    bams="${OUT_bwa_fna}/${line}"*_mapped_sort_NOdup_sort.bam
    echo "samtools merge "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged.bam" "${bams}"" >> $bash_out
    # index
    echo "samtools index "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged.bam"" >> $bash_out
    # sort
    echo "samtools sort "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged.bam" -o "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged_sort.bam"" >> $bash_out
    # index
    echo "samtools index "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged_sort.bam"" >> $bash_out
    # make a sam file
    echo "samtools view -h -o "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged_sort.sam" "${OUT_bwa_fna}/${line}_mapped_sort_NOdup_sort_merged_sort.bam"" >> $bash_out

    # move annotated file
    #OUT_gbk_annotated="${OUT_gbk}/output/evidence/annotated.gd"
    #essentials_gbk_annotated="${OUT_gbk_essentials}/annotated.gd"
    #echo "cp ${OUT_gbk_annotated} ${essentials_gbk_annotated}" >> $bash_out

    # move evidence file
    #OUT_gbk_evidence="${OUT_gbk}/output/evidence/evidence.gd"
    #essentials_gbk_evidence="${OUT_gbk_essentials}/evidence.gd"
    #echo "cp ${OUT_gbk_evidence} ${essentials_gbk_evidence}" >> $bash_out

    # move output file
    #OUT_gbk_output="${OUT_gbk}/output/output.gd"
    #essentials_gbk_output="${OUT_gbk_essentials}/output.gd"
    #echo "cp ${OUT_gbk_output} ${essentials_gbk_output}" >> $bash_out


    # get coverage from bam
    #bam="${OUT}/data/reference.bam"
    #coverage="${OUT}/data/coverage.txt"
    #echo "samtools mpileup ${bam} > ${coverage}" >> $bash_out

    # calculate dist of coverage
    # this should be done as a seperate command
    #coverage_clean="${OUT_essentials}/coverage_clean.txt"
    #mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}"
    #echo "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}"
    #OUT_coverage="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}/coverage_clean.txt"
    #echo "python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -c -i ${coverage} -o ${OUT_coverage}" >> $bash_out
    qsub $bash_out
  fi
done
