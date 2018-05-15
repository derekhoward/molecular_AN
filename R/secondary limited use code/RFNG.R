#short bit of code to get RFNG rankings
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)

geneOfInterest <- "RFNG"
rfng <- read_tsv("./data/processed/adult_brainarea_vs_genes_exp_reannotator.tsv")
z <- melt(rfng) %>% as_tibble()
z %<>% group_by(variable) %>% mutate(rank = rank(value))

#z %>% filter(variable == "lateral parabrachial nucleus")

z <- filter(z, gene_symbol==geneOfInterest)
#z <- melt(z) %>% as_tibble()
z %>% arrange(-rank)
threshold <- z %>% filter(variable == 'lateral parabrachial nucleus') %>% .$value

rfng <- read_tsv("./data/processed/adult_brainarea_vs_genes_exp_reannotator.tsv")

rfng %>% select(gene_symbol, `lateral parabrachial nucleus`) %>% arrange(-`lateral parabrachial nucleus`) 

rfng %>% select(gene_symbol, `lateral parabrachial nucleus`) %>% filter(`lateral parabrachial nucleus` > threshold)

