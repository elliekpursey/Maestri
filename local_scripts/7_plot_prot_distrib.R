library(tidyverse)
library(MetBrewer)
library(ggpubr)

# set working dir
root_directory <- "" # where your results files are

setwd(root_directory)

best_solutions_df <- read_csv('best_solution.csv', col_names=c('replicon', 'hit_id', 'gene_name', 'hit_pos',
                                                      'model_fqn', 'sys_id',	'sys_loci',	'locus_num',
                                                      'sys_wholeness',	'sys_score',	'sys_occ',	'hit_gene_ref',
                                                      'hit_status',	'hit_seq_len',	'hit_i_eval',	'hit_score',
                                                      'hit_profile_cov',	'hit_seq_cov', 'hit_begin_match',	'hit_end_match',
                                                      'used_in', 'file'))

species_df <- read_csv("species_accessions_refseq_bacteria.csv") %>%
  rename(genome = accession) %>%
  separate(col=species, into=c("genus", "species"), sep=" ")

tidy_df <- best_solutions_df %>%
  separate(col=file, into=c("results", "macsy_results", "genome"), sep="/") %>%
  mutate(full_genome = genome) %>%
  separate(col=genome, into=c("GCF", "id_no"), sep="_") %>%
  unite(genome, GCF:id_no, sep = '_', remove = FALSE) %>%
  dplyr::select('replicon', 'gene_name', 'hit_id', 'hit_status', 'hit_begin_match', 'hit_end_match', 'genome', 'full_genome') %>%
  full_join(., species_df) %>%
  mutate(gene_name = str_replace(gene_name, "prot_", ""))

# presence of each gene in systems

prot.colours <- c("1" = "#DF9ED4", "2" = "#D16284", "3" = "#D16258", "4" = "#E5BA60",
                                "5" = "#98B46C", "6" = "#44927A", "7" = "#3E5F90", "8" = "#5C4699")

count_genes <- tidy_df %>%
  select(full_genome, gene_name) %>%
  unique() %>%
  count(gene_name) %>%
  mutate(total_genomes = 422) %>% # always check this is the total number you have
  mutate(percent = (n/total_genomes) * 100)

summary_plot <- ggplot(count_genes, aes(x = gene_name, y = percent, fill=gene_name)) +
  theme_light() +
  theme(legend.position = "False") +
  geom_col(position = "dodge") +
  labs(x="Protein in Maestri operon", y="% of systems where protein detected") +
  scale_fill_manual(values=prot.colours)

ggsave("gene_distrib_plot.jpg", summary_plot, width=18, height=10, dpi=300, units="cm")

# gene distributions for selected genera

plot_protein_distrib <- function(genus_name){
  
  genus_df <- tidy_df %>%
    dplyr::filter(genus == genus_name)
  
  total_genomes_length <- length(unique(genus_df$genome)) 
  
  count_genes_genus <- genus_df %>%
    select(full_genome, gene_name) %>%
    unique() %>%
    count(gene_name) %>%
    mutate(total_genomes = total_genomes_length) %>% 
    mutate(percent = (n/total_genomes) * 100) %>%
    mutate(gene_name = str_replace(gene_name, "prot_", ""))
  
  ylab <- bquote('% in '~italic(.(genus_name))) 
                 
  plot <- ggplot(count_genes_genus, aes(x = gene_name, y = percent, fill=gene_name)) +
    theme_light() +
    theme(legend.position = "False") +
    geom_col(position = "dodge") +
    labs(x="Protein in Maestri operon", y=ylab) +
    scale_fill_manual(values=prot.colours) 
  
  assign(paste0(genus_name, "_plot", sep=""), plot,  envir = parent.frame())
  
}
  
genus_list <- c("Escherichia", "Pseudomonas", "Ralstonia", "Streptomyces", "Klebsiella", "Vibrio")

for(genus_name in genus_list){
  plot_protein_distrib(genus_name)
}

prot_distrib_arrange <- ggarrange(Escherichia_plot, Pseudomonas_plot, Ralstonia_plot, Streptomyces_plot, Klebsiella_plot, Vibrio_plot,
          nrow = 3, 
          ncol = 2,
          labels = c("B", "C", "D", "E", "F", "G"), 
          font.label = list(size = 14))

final_plot <- ggarrange(summary_plot, prot_distrib_arrange,
                        nrow=1,
                        ncol=2,
                        widths = c(1.5,1),
                        labels = c("A", ""),
                        font.label = list(size = 14))

ggsave("example_and_all_distrib_plots.jpg", final_plot, width=35, height=20, dpi=300, units="cm")
ggsave("example_and_all_distrib_plots.svg", final_plot, width=35, height=20, dpi=300, units="cm")

