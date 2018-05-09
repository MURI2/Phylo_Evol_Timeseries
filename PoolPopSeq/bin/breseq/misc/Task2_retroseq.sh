echo "#!/bin/bash" > Task3${i}_EcoliRetroSeq.sh
echo "" >> Task3${i}_EcoliRetroSeq.sh
echo "#$ -l vmem=50gb walltime=2:00:00 " >> Task3${i}_EcoliRetroSeq.sh
echo "" >> Task3${i}_EcoliRetroSeq.sh
echo "cd /N/dc2/projects/muri2/Task3/${i}/"  >>Task3${i}_EcoliRetroSeq.sh
echo "mkdir E_coli_InsSeq" >> Task3${i}_EcoliRetroSeq.sh
echo "cd E_coli_InsSeq" >>Task3${i}_EcoliRetroSeq.sh
echo "time bwa mem /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna ../Sample_${i}_R1_trimmed.fastq ../Sample_${i}_R2_trimmed.fastq > Sample_${i}.sam" >> Task3${i}_EcoliRetroSeq.sh
echo "samtools view -bS -T /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna Sample_${i}.sam > Sample_${i}.bam" >> Task3${i}_EcoliRetroSeq.sh
echo "samtools sort Sample_${i}.bam Sample_${i}.sorted" >> Task3${i}_EcoliRetroSeq.sh
echo "samtools index Sample_${i}.sorted.bam" >> Task3${i}_EcoliRetroSeq.sh
echo "perl /N/dc2/projects/muri2/Tools/RetroSeq/bin/retroseq.pl -discover -eref /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/Ecoli_IS_RefFile.txt -bam Sample_${i}.sorted.bam -output Sample_${i}.IS.Reads -align" >> Task3${i}_EcoliRetroSeq.sh
echo "perl /N/dc2/projects/muri2/Tools/RetroSeq/bin/retroseq.pl -call -bam Sample_${i}.sorted.bam -ref /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna -output Sample_${i}.Ecoli.IS -input Sample_${i}.IS.Reads -hets" >>Task3${i}_EcoliRetroSeq.sh
echo "" >> Task3${i}_EcoliRetroSeq.sh
echo "exit" >> Task3${i}_EcoliRetroSeq.sh
