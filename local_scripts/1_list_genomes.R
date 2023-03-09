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

tidy_df <- best_solutions_df %>%
  separate(col=file, into=c("results", "macsy_results", "genome"), sep="/") %>%
  mutate(full_genome = genome) %>%
  separate(col=genome, into=c("GCF", "id_no"), sep="_") %>%
  unite(genome, GCF:id_no, sep = '_', remove = FALSE) %>%
  dplyr::select('replicon', 'gene_name', 'hit_id', 'hit_status', 'hit_begin_match', 'hit_end_match', 'genome', 'full_genome')

write_csv(tidy_df, file = "tidy_systems_df.csv")

genome_names <- tidy_df %>%
  separate(col=genome, into=c("GCF", "id_no"), sep="_") %>%
  unite(genome, GCF:id_no, sep = '_', remove = FALSE) %>%
  dplyr::select(genome) %>%
  unique()

write_csv(genome_names, file = "genome_list.csv")
write_delim(genome_names, "genome_list.txt", delim="\n", col_names=FALSE)

# count number of proteins per system
count_files <- tidy_df %>%
  count(genome, name="number_of_proteins_in_system")

count_system_gene_totals <- count_files %>%
  count(number_of_proteins_in_system, name="number_of_systems") %>%
  mutate(number_of_proteins_in_system = as.factor(number_of_proteins_in_system))

number_of_proteins_per_system <- ggplot(count_system_gene_totals, aes(x=number_of_proteins_in_system, 
                                                                      y=number_of_systems)) +
  geom_col() +
  theme_light() +
  labs(x="Number of proteins detected in Maestri", y="Number of systems") 

ggsave("Number_of_proteins_per_system.svg", plot=number_of_proteins_per_system, width=10, height=10, units="cm", dpi=300)
  
