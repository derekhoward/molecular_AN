library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)
library(homologene)

batchesToUse <- c("B1", "B2") #either B1 or B2 - batches with hunger and normal cells

#batchesToUse <- c("B0", "B1", "B2", "B3", "B4") #either B1 or B2 - batches with hunger and normal cells

#can break Rstudio when it trys to show the matrix in the environment panel, run outside of rstudio to process it
load("./data/raw/single_cell/GSE87544_1443737Cells.Expresssion.Matrix.log_tpm+1_.renamed.RData", verbose = T)
GSE87544 <- Expresssion_Matrix_unfiltered
Expresssion_Matrix_unfiltered <- NULL
dim(GSE87544)
#GSE87544 <- head(GSE87544,n=1000) #for testing 

targetGeneList <- "allInHomologene"
allHuman <- unique(tbl_df(homologeneData) %>% filter(Taxonomy == 9606) %>% .$Gene.Symbol)
targetSymbols <- unique(human2mouse(allHuman)$mouseGene)

#cut down the matrix to genes of interest so easier to handle
GSE87544 <- GSE87544[intersect(rownames(GSE87544),targetSymbols),]

#some duplicated code with the cell-type specific code
cellClusters <- read_csv("./data/raw/single_cell/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz")
cellClusters %<>% tbl_df() %>% rename(cell_id = X1)
cellClusters %<>% mutate(batch = substr(cell_id, 1, 2), group = substr(cell_id, 16, length(cell_id)))
cellClusters %<>% mutate(group_combined = gsub("[12.]","", group))
cellClusters %<>% select(-Cell_ID_temp)

cellClusters %>% group_by(batch) %>% summarize(n=n())
cellClusters %>% group_by(batch, group_combined) %>% summarize(n=n())
cellClusters %>% group_by(batch, group) %>% summarize(n=n())
cellClusters %>% group_by(group) %>% summarize(n=n())
cellClusters %>% group_by(group_combined) %>% summarize(n=n())

#filter batches that have only one condition
cellClusters %<>% filter(batch %in% batchesToUse)
cellClusters %>% group_by(batch, group_combined) %>% summarize(n=n())

clusterSizes <- cellClusters %>% group_by(SVM_clusterID) %>% summarize(cells=n()) %>% arrange(desc(cells))

GSE87544$geneSymbol <- rownames(GSE87544)

GSE87544Melted <- reshape2::melt(GSE87544,factorsAsStrings = TRUE, id.vars=c("geneSymbol"), variable.name = "cell_id", value.name = "log1Expression")

GSE87544Melted <- tbl_df(GSE87544Melted)
GSE87544Melted

#join with cell information
GSE87544Melted <- inner_join(select(cellClusters, -group), GSE87544Melted)
#write_csv(GSE87544Melted, "./results/cell types/Rps26.GSE87544.csv") #for writing single gene

#filter all zero genes
zeroGenes <- GSE87544Melted %>% group_by(geneSymbol, batch) %>% summarize(maxExp = max(log1Expression)) %>% filter(maxExp == 0)
GSE87544Melted %<>% filter(!geneSymbol %in% zeroGenes$geneSymbol)

#find cell types with no hunger/normal cells
limitedGroups <- GSE87544Melted %>% group_by(geneSymbol, SVM_clusterID, batch) %>% summarize(
  n = n(), 
  groups = length(unique(group_combined))
) %>% filter(groups < 2) %>% ungroup() %>% select(SVM_clusterID, batch)
unique(limitedGroups$SVM_clusterID)

#split by cell type
hungerPerGeneByType <- GSE87544Melted %>% group_by(geneSymbol, SVM_clusterID, batch) %>% anti_join(limitedGroups) %>% summarize(
  n = n(),
  nonzero = sum(log1Expression > 0 ), 
  pValue = wilcox.test(log1Expression ~ group_combined)$p.value, 
  medianExp = median(log1Expression)
) %>% arrange(pValue) %>% ungroup() 
hungerPerGeneByType %<>% filter(!is.na(pValue))
hungerPerGeneByType %<>% group_by(batch) %>% mutate(pValue.adj = p.adjust(pValue, method="fdr")) %>% print()
filter(hungerPerGeneByType, pValue.adj < 0.05)
strongestCellType <- hungerPerGeneByType %>% arrange(pValue) %>% group_by(geneSymbol, batch) %>% summarize(strongestCellType = first(SVM_clusterID))

meansPerGroupAndType <- GSE87544Melted %>% group_by(geneSymbol, SVM_clusterID, group_combined, batch) %>% summarize(
  meanExp = mean(log1Expression),
)
meansPerGroupAndType %<>% spread( group_combined, meanExp)
hungerPerGeneByType <- inner_join(meansPerGroupAndType, hungerPerGeneByType)


write_csv(hungerPerGeneByType, paste0("./results/cell types/", targetGeneList, ".GSE87544.hungerPerGene.wilcox.byType.csv"))

#don't split by cell type
hungerPerGene <- GSE87544Melted %>% group_by(geneSymbol, batch) %>% summarize(
  n = n(),
  nonzero = sum(log1Expression > 0 ), 
  pValue = wilcox.test(log1Expression ~ group_combined)$p.value, 
  medianExp = median(log1Expression)
) %>% arrange(pValue) %>% group_by(batch) %>% mutate(pValue.adj = p.adjust(pValue, method="fdr"))

#add in strongest cell type
hungerPerGene <- inner_join(strongestCellType, hungerPerGene)

mediansPerGroup <- GSE87544Melted %>% group_by(geneSymbol, group_combined, batch) %>% summarize(
  n=n(),
  meanExp = mean(log1Expression),
  nonzero = sum(log1Expression > 0)
)
mediansPerGroup <- select(mediansPerGroup,-nonzero, -n) %>% spread( group_combined, meanExp)
hungerPerGene <- inner_join(mediansPerGroup, hungerPerGene)

hungerPerGene %<>% mutate(direction = sign(Hungry - Normal)) %>% arrange(pValue)  %>% print()
filter(hungerPerGene, pValue.adj < 0.05) %>% print() %>% group_by(geneSymbol) %>% summarize(n=n(), directionSum = sum(direction)) %>% arrange(desc(n))

write_csv(hungerPerGene, paste0("./results/cell types/", targetGeneList, ".GSE87544.hungerPerGene.wilcox.csv"))

#for visualizing a single gene
#forPlot <- filter(GSE87544Melted, geneSymbol == "Pfn1")
#ggplot(forPlot, aes(y = log1Expression, x = SVM_clusterID, color = group_combined)) + geom_violin()

#ggplot(forPlot, aes(y = log1Expression, x = group_combined)) + geom_violin()
