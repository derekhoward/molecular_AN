library(broom)
library(tidyr)
library(magrittr)
library(homologene)
library(readr)
library(dplyr) 
library(reshape2)

batchesToUse <- c("b5", "b6") #either b5 or b6 - batches with hunger and normal cells

GSETableHeader <- read_tsv("./data/raw/single_cell/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt.gz", n_max=1, col_names = F)
GSETable <- read_tsv("./data/raw/single_cell/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt.gz", col_names = F, skip = 1)
colnames(GSETable) <- c("geneSymbol", GSETableHeader[1,])
#GSETable %<>% filter(geneSymbol == "Rps26") #for single cell data filtering
#GSETable <- head(GSETable, n=1000) #for testing only
dim(GSETable)

targetGeneList <- "allInHomologene"
allHuman <- unique(tbl_df(homologeneData) %>% filter(Taxonomy == 9606) %>% .$Gene.Symbol)
targetSymbols <- unique(human2mouse(allHuman)$mouseGene)

GSETableTargetted <- GSETable %>% filter(geneSymbol %in% targetSymbols)
GSETable <- NULL
#some replicated code
#get cell info
cellTable <- read_tsv("./data/raw/single_cell/GSE93374_cell_metadata.txt.gz")
colnames(cellTable) <- gsub("[0-9]+[.]","", colnames(cellTable))

##ch10 is 10% calories from fat
cellTable %>% group_by(batches, FvF) %>% summarize(n=n())

#cell cluster lookup table to handle Fibroblasts to VLMCs changes
lookup <- read_csv("./data/raw/single_cell/GSE93374_cell labels from supplement spreadsheets.txt", col_names = F) %>% rename(subcluster = X1)
lookup %<>% mutate(clust_all_micro =gsub("[.].*", "", subcluster), clusterName = gsub(".*[.]", "", subcluster))

cellTable %<>% select(ID, group, sex, Sex_pred, FvF, clust_all_micro, batches)

cellTable <- inner_join(cellTable, select(lookup, clust_all_micro, clusterName)) %>% select(-clust_all_micro)

cellTable %<>% filter(batches %in% batchesToUse)
cellTable %<>% filter( sex == "F" | Sex_pred == "F") #chen was all females, known sex is for only batch 6 of Campbell


GSETableTargetted <- melt(GSETableTargetted, factorsAsStrings = TRUE, id.vars=c("geneSymbol"), variable.name = "ID", value.name = "log1Expression")
GSETableTargetted <- tbl_df(GSETableTargetted)
dim(GSETableTargetted)

GSETableTargetted <- inner_join(GSETableTargetted, cellTable)
#write_csv(GSETableTargetted, "./results/cell types/Rps26.GSE93374.csv")

GSETableTargetted %>% group_by(batches, sex) %>% summarize(n=length(unique(ID)))
GSETableTargetted %>% group_by(batches, Sex_pred) %>% summarize(n=length(unique(ID)))


noVarianceGenes <- GSETableTargetted %>% group_by(geneSymbol, batches) %>% summarize(maxExp = max(log1Expression), minExp = min(log1Expression)) %>% filter(maxExp - minExp == 0)
GSETableTargetted %<>% filter(!geneSymbol %in% noVarianceGenes$geneSymbol)

length(unique(GSETableTargetted$geneSymbol))

#ch10 is 10% calories from fat
GSETableTargetted %>% group_by(batches, FvF) %>% summarize(n=length(unique(ID)))

hungerPerGeneGSE93374 <- GSETableTargetted %>% group_by(geneSymbol, batches) %>% summarize(
  n = n(),
  nonzero = sum(log1Expression > 0 ), 
  pValue = wilcox.test(log1Expression ~ FvF)$p.value, 
  medianExp = median(log1Expression)
) %>% arrange(pValue)

meansPerGroup <- GSETableTargetted %>% group_by(geneSymbol, FvF, batches) %>% summarize(
  n=n(),
  meanExp = mean(log1Expression),
  nonzero = sum(log1Expression > 0)
)
meansPerGroup <- select(meansPerGroup,-nonzero, -n) %>% spread( FvF, meanExp)
hungerPerGeneGSE93374 <- inner_join(meansPerGroup, hungerPerGeneGSE93374)

hungerPerGeneGSE93374 %<>% group_by(batches) %>% mutate(pValue.adj = p.adjust(pValue, method="fdr")) %>% arrange(pValue)

hungerPerGeneGSE93374 %<>% mutate(direction = if_else(batches == "b5", sign(Fast - Refed), sign(Fast - Fed)))

nrow(hungerPerGeneGSE93374 %>% filter(pValue < 0.05))

write_csv(hungerPerGeneGSE93374, paste0("./results/cell types/", targetGeneList, ".GSE93374.hungerPerGene.wilcox.csv"))
