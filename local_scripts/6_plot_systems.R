library(tidyverse)
library(gggenes)
library(ghibli)
library(svglite)

# set working dir
root_directory <- "" # where your results files are

setwd(root_directory)


# assign colours to genes in Maestri

# diagram of Maestri locus

Maestri_system <- as.data.frame(list(genome_locus = c('Maestri', 'Maestri', 'Maestri', 'Maestri', 'Maestri', 'Maestri', 'Maestri', 'Maestri'),
                                   strand = c(-1, -1, -1, -1, -1, -1, -1, -1),
                                   gene_hit = c("8", "7", "6", "5", "4", "3", "2", "1"),
                                   hit_start_loc = c(-135866, -141457, -143130, -147307, -148716, -149369, -150703, -152738),
                                   hit_end_loc = c(-141460, -143130, -147305, -148719, -149372, -150706, -152745, -152953),
                                   genus = c('Pseudomonas', 'Pseudomonas', 'Pseudomonas', 'Pseudomonas', 'Pseudomonas', 'Pseudomonas', 'Pseudomonas', 'Pseudomonas'),
                                   species = c('aeruginosa', 'aeruginosa', 'aeruginosa', 'aeruginosa', 'aeruginosa', 'aeruginosa', 'aeruginosa', 'aeruginosa')))

                                   

# plot maestri system

plot_maestri<- ggplot(Maestri_system, aes(xmin = hit_start_loc, xmax = hit_end_loc, y=genome_locus, label=gene_hit, forward = strand)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), alpha = 1, colour="black") +
  geom_gene_label() + 
  theme_genes() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position="none")

ggsave("maestri_locus_diagram.svg", plot=plot_maestri, units="cm", width=25, height=5, dpi=300)

# plot maestri system with domain hits

HHpred_hits <- read_csv('hhpred_domains.csv') %>%
  mutate(gene_hit = as.factor(gene_hit)) %>%
  full_join(., Maestri_system) %>%
  mutate(subgene_min = hit_start_loc - (domain_begin*3)) %>% # multiply by three to convert aa length to nucleotide length
  mutate(subgene_max = hit_start_loc - (domain_end*3)) 

HHpred_hits$domain_type <- factor(HHpred_hits$domain_type, levels = c("DNA-binding, HTH, regulator", "N6 methyltransferase",
                                                                      "ATPase", "DNA-binding domain", "DUF4276", "type I R-M specificity subunit",
                                                                      "protein kinase", "PT-dependent restriction protein", "helicase, translocase"))

plot_maestri_domains <- ggplot(HHpred_hits, aes(xmin = hit_start_loc, xmax = hit_end_loc, y=genome_locus, label=domain_type, fill = domain_type, forward = strand)) +
  geom_gene_label() + 
  theme_genes() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), alpha = 1, fill="white", colour="black") +
  geom_subgene_arrow(data = HHpred_hits,
                   aes(xmin = hit_start_loc, xmax = hit_end_loc, y=genome_locus, fill = domain_type,
                       xsubmin = subgene_min, xsubmax = subgene_max), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), color="black") +
  scale_fill_brewer(palette ="Set3", name = "Protein and domain hits") 

ggsave("maestri_domain_hits_diagram.svg", plot=plot_maestri_domains, units="cm", width=30, height=10, dpi=300)


