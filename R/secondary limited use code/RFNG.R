geneOfInterest <- "RFNG"
rfng <- read_tsv("./data/processed/adult_brainarea_vs_genes_exp_reannotator.tsv")
z <- filter(rfng, gene_symbol==geneOfInterest)
z <- melt(z) %>% as_tibble()
z %>% arrange(-value)
threshold <- z %>% filter(variable == 'lateral parabrachial nucleus') %>% .$value

rfng %>% select(gene_symbol, `lateral parabrachial nucleus`) %>% arrange(-`lateral parabrachial nucleus`) 

rfng %>% select(gene_symbol, `lateral parabrachial nucleus`) %>% filter(`lateral parabrachial nucleus` > threshold)

