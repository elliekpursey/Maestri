from Bio import SeqIO
import os 
import glob
import gzip
import shutil

file = snakemake.input[0]
genome = snakemake.wildcards['genome']

# gunzip gbff file
with gzip.open(file, 'rb') as f_in:
    gunzip_file = "resources/all_genomes/" + genome + ".gbff"
    with open(gunzip_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# initialise empty list for CDS entries
all_entries = []

# add all CDS to list and write to file
with open(gunzip_file, 'r') as gbk:
    protein_seqs = SeqIO.InsdcIO.GenBankCdsFeatureIterator(gbk)
    for cds in protein_seqs:
        if cds.seq is not None:
            cds.id = cds.name
            cds.description = ''
            all_entries.append(cds)

out_fasta = "results/protein_fastas/" + genome + ".fasta"
SeqIO.write(all_entries, out_fasta, 'fasta')     

# gzip gbff file and remove uncompressed version
with open(gunzip_file, 'rb') as f_in, gzip.open(file, 'wb') as f_out:
    f_out.writelines(f_in)

os.remove(gunzip_file)

# gzip fasta file and remove uncompressed version
with open(out_fasta, 'rb') as f_in, gzip.open(snakemake.output[0], 'wb') as f_out:
    f_out.writelines(f_in)

os.remove(out_fasta)



