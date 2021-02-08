[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4517573.svg)](https://doi.org/10.5281/zenodo.4517573)



# Phylo_Evol_Timeseries


Repository for code associated with the preprint: "Molecular evolutionary dynamics of energy limited microorganisms"



## Dependencies
All code was written in Python 3.6. The conda environment is detailed in `environment.yml`. Set up the repo under a folder named Github: `~/GitHub/Phylo_Evol_Timeseries/`.

## Getting the data

Assembled genomes are available on NCBI and accession numbers listed in the preprint. Raw reads are available on SRA under BioProject PRJNA639414. All other data is on Zenodo data repository doi:10.5281/zenodo.4517573.


## Running the analyses

Put all reads from SRA under the folder `Phylo_Evol_Timeseries/data/illumina_runs`.




Just run `run_everything.sh`. This script contains commands for all data processing and analyses, including mutation calling. You will need to rework the paths in the script to get it working on your machine. All bash scripts are written to be run on a cluster. I do not recommend running them on your local machine.
