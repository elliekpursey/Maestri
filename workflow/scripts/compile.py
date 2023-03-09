#!/usr/bin/env python3

import pandas as pd
import glob

input_list = glob.glob("results/macsy_results/*")

# best solution
for folder in input_list:
    try:
        df = pd.read_table(folder + "/best_solution.tsv", sep="\t", skiprows=3)
        df['genome'] = folder
        df.to_csv(snakemake.output[0], mode='a', header=False, index=False) # append to output csv
    except pd.errors.EmptyDataError:
        print("No data")

# all best solutions
for folder in input_list:
    try:
        df = pd.read_table(folder + "/all_best_solutions.tsv", sep="\t", skiprows=3)
        df['genome'] = folder
        df.to_csv(snakemake.output[1], mode='a', header=False, index=False) # append to output csv
    except pd.errors.EmptyDataError:
        print("No data")

# all systems
for folder in input_list:
    try:
        df = pd.read_table(folder + "/all_systems.tsv", sep="\t", skiprows=3)
        df['genome'] = folder
        df.to_csv(snakemake.output[2], mode='a', header=False, index=False) # append to output csv
    except pd.errors.EmptyDataError:
        print("No data")
        