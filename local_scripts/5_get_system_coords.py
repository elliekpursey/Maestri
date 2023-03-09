#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import pandas as pd

results_dir = "" # where your results files are

os.chdir(results_dir)

df = pd.read_csv("tidy_systems_df.csv")

genome_grouped = df.groupby('full_genome')

model_hits = []
loci = []
genomes = []
proteins = []
locations = []
products = []
hit_begins = []
hit_ends = []

for name, group in genome_grouped:
    try:
        for record in SeqIO.parse("genomes_with_systems/" + name + '.gbff', "genbank"):
            locus = record.name
            for feature in record.features:
                if (feature.type == "CDS") & ("locus_tag" in feature.qualifiers):
                    product = str(feature.qualifiers['product'][0])
                    protein_id = str(feature.qualifiers['locus_tag'][0])
                    loc = feature.location
                    for row_index, row in group.iterrows():
                        if row['hit_id'] == protein_id: 
                            print('adding ' + protein_id)
                            hit_begin = row['hit_begin_match']
                            hit_end = row['hit_end_match']
                            model_hit = row['gene_name']
                            model_hits.append(model_hit)
                            loci.append(locus)
                            genomes.append(name)
                            proteins.append(protein_id)
                            products.append(product)
                            locations.append(loc)
                            hit_begins.append(hit_begin)
                            hit_ends.append(hit_end)
    except Exception as e:
        print(e)
        
df_final = pd.DataFrame(list(zip(model_hits, loci, genomes, proteins,locations,products, hit_begins, hit_ends)),
               columns =['gene_hit','locus','genome', 'protein', 'location', 'product', 'hit_begin', 'hit_end'])

df_final.to_csv('system_coords.csv', index=False)
