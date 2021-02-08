#!/bin/bash



# clean raw reads
python Python/clean_file_names.py

bash bash/rename_reads.sh

bash bash/reads_clean.sh

bash bash/breseq.sh

bash bash/breseq_jc.sh

bash bash/rebreseq.sh

bash bash/create_timecourse.sh

bash bash/create_merged_timecourse.sh

bash bash/annotate_pvalues.sh




# process mutation trajectories
bash raxml.sh

python process_cluster_output.py

python calculate_convergence_matrices.py

python calculate_parallel_genes_table.py



# analyze mutation trajectories

python plot_allele_frequencies.py

python plot_convergence_genes.py

python plot_diversity.py

python analyze_mutation_specrta.py

python analyze_extinctions.py

python plot_dnds.py

python plot_fmax_parallelism.py

python plot_genome_wide_parallelism_and_divergence.py

python plot_multiplicity.py

python plot_parallel_genes.py

python plot_parallel_pathways_and_get_table.py

python plot_per_taxon_divergence.py

python plot_rate.py
