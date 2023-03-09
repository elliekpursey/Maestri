library(tidyverse)
library(MetBrewer)

# set working dir
root_directory <- "" # where your results files are

setwd(root_directory)

data_raw <- read_csv('best_solution.csv', col_names=c('replicon', 'hit_id', 'gene_name', 'hit_pos',
                                                      'model_fqn', 'sys_id',	'sys_loci',	'locus_num',
                                                      'sys_wholeness',	'sys_score',	'sys_occ',	'hit_gene_ref',
                                                      'hit_status',	'hit_seq_len',	'hit_i_eval',	'hit_score',
                                                      'hit_profile_cov',	'hit_seq_cov', 'hit_begin_match',	'hit_end_match',
                                                      'used_in', 'file'))

data_tidy <- data_raw %>%
  separate(col=file, into=c("results", "macsy_results", "genome"), sep="/") %>%
  select('replicon', 'gene_name', 'hit_id', 'hit_status', 'hit_begin_match', 'hit_end_match', 'genome')

count_genes <- data_tidy %>%
  select(genome, gene_name) %>%
  unique() %>%
  count(gene_name) %>%
  mutate(total_genomes = 172366) %>% # always check this is the total number you have
  mutate(percent = (n/total_genomes) * 100) %>%
  mutate(gene_name = str_replace(gene_name, "prot_", ""))

summary_plot <- ggplot(count_genes, aes(x = gene_name, y = percent, fill=gene_name)) +
  theme_light() +
  theme(legend.position = "False") +
  geom_col(position = "dodge") +
  labs(x="Gene in Maestri operon", y="% of genomes where gene detected") +
  scale_fill_manual(values=met.brewer("Klimt", 9), name = "Number of genes in model") 

ggsave("gene_distrib_plot.jpg", summary_plot, width=18, height=10, units="cm")

