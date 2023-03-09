library(tidyverse)

# set working dir
root_directory <- "" # where your results files are

setwd(root_directory)

# read in df of all solutions

best_solutions_df <- read_csv('best_solution.csv', col_names=c('replicon', 'hit_id', 'gene_name', 'hit_pos', 
                                                                        'model_fqn', 'sys_id',	'sys_loci',	'locus_num',	
                                                                        'sys_wholeness',	'sys_score',	'sys_occ',	'hit_gene_ref',
                                                                        'hit_status',	'hit_seq_len',	'hit_i_eval',	'hit_score',
                                                                        'hit_profile_cov',	'hit_seq_cov', 'hit_begin_match',	'hit_end_match',	
                                                                        'used_in', 'file'))

tidy_best_solutions_df <- best_solutions_df %>%
  separate(col=file, into=c("results", "macsy_results", "genome"), sep="/") %>%
  separate(col=genome, into=c("GCF", "id_no"), sep="_") %>%
  unite(genome, GCF:id_no, sep = '_', remove = FALSE) %>%
  dplyr::select(genome, gene_name, sys_wholeness, sys_score)

# read in species df 

species_df <- read_csv("species_accessions_refseq_bacteria.csv") %>%
  rename(genome = accession)

unique_species <- species_df %>%
  dplyr::select(species) %>%
  unique()

# join dataframes

species_systems_df <- tidy_best_solutions_df %>%
  full_join(., species_df) 

# count species and plot

count_plot_genus_df <- species_systems_df %>%
  dplyr::select(species, genome) %>%
  unique() %>%
  separate(col=species, into=c("genus", "species"), sep=" ") %>%
  count(genus) %>%
  arrange(-n) %>%                              
  mutate(genus = factor(genus, genus))

genus_plot <- ggplot(count_plot_genus_df, aes(x = genus, y = n)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_col(position = "dodge") +
  labs(x="Genus", y="No. of systems") 

ggsave("systems_by_genus.svg", plot=genus_plot, width=18, height=15, dpi=300, units="cm")

