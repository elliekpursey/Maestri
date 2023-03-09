#!/bin/bash 

cd # where your results files are

# download final_model_genomes

mkdir genomes_with_systems
cd genomes_with_systems

ncbi-genome-download bacteria -F 'genbank' --flat-output -A ../genome_list.txt --verbose

for file in *.gz; do gunzip $file; done
