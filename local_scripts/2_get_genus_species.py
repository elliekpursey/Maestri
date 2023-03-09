#!/usr/bin/env python3

from Bio import Entrez
from Bio import SeqIO
import pandas as pd

# results dir
results_dir = ' ' # location of your results

# read in list of genomes
genomes = pd.read_csv(results_dir + "/genome_list.csv", header=0)
genome_list = list(genomes.genome)

Entrez.email = ' ' # your email here
species_list = []
species_taxid_list = []
id_list = []

for id in genome_list:
    print(id)
    handle = Entrez.esearch(db="assembly", term=id, retmax='1')
    record = Entrez.read(handle, validate=False)
    esummary_handle = Entrez.esummary(db="assembly", id=record['IdList'], report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    species = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
    species_taxid = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesTaxid']
    species_list.append(species)
    species_taxid_list.append(species_taxid)
    id_list.append(id)

df = pd.DataFrame(list(zip(species_list, species_taxid_list, id_list)),
            columns =['species','species_taxid','accession'])

df.to_csv(results_dir + '/species_accessions_refseq_bacteria.csv', index=False)

