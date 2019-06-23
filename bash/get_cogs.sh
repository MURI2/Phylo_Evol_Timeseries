

conda create -n anvio5_test -c bioconda -c conda-forge anvio=5.5.0 python=3.6

#for d in /Users/WRShoemaker/GitHub/LTDE/data/genomes/nanopore_hybrid_annotated/*; do
#    strain="$(echo "$d" | cut -d "/" -f9-9 )"
#    echo $strain
#    REF_OUT="~/LTDE/data/genomes/nanopore_hybrid_annotated_cogs/${strain}_reformat"
#    anvi-script-reformat-fasta $d/FCE86-Genome.fna -o $REF_OUT.fna -l 0 --simplify-names
#    anvi-gen-contigs-database -f $REF_OUT.fna -o $REF_OUT.db
#    anvi-run-ncbi-cogs -c $REF_OUT.db --num-threads 20
#    anvi-export-functions -c $REF_OUT.db -o $REF_OUT.txt
#done
