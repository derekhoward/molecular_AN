library(reshape2)
library(readr)
library(homologene)
library(metap)
library(dplyr)
library(magrittr)
library(tidyr)

targetGeneLists <- c("Negraes et al. Table S5.hypenFixed.mouse.txt", "Duncan et al. rs4622308.hypenFixed.mouse.txt", "Lutter et al. Table S3.RestrictedEating.hypenFixed.mouse.txt", "Lutter et al. Table S4.BingeEating.hypenFixed.mouse.txt")

GSE87544_Chen <- read_csv(paste0("./results/cell types/allInHomologene.GSE87544.hungerPerGene.wilcox.csv"))
GSE93374_Cambpell <- read_csv(paste0("./results/cell types/allInHomologene.GSE93374.hungerPerGene.wilcox.csv"))
GSE93374_Cambpell %<>% rename(Hungry = Fast, Normal = Fed, batch = batches)
GSE93374_Cambpell %<>% filter(batch != "b5") %>% select(-Refed)

bothDataSets <- bind_rows(GSE87544_Chen, GSE93374_Cambpell) 
bothDataSets %<>% mutate(pvalue.pos = if_else(sign(direction) > 0, pValue, 1 - pValue )) 
bothDataSets %<>% mutate(pvalue.neg = if_else(sign(direction) < 0, pValue, 1 - pValue ))

bothDataSetsSummary <- bothDataSets %>% group_by(geneSymbol) %>% filter(n() ==3 ) %>% summarize(
  directionSum = sum(direction), 
  n = n(), 
  metaP.pos = sumlog(pvalue.pos)$p, 
  metaP.neg = sumlog(pvalue.neg)$p
)
bothDataSetsSummary %<>% mutate(metaP.pos.adj = p.adjust(metaP.pos, method="fdr"))
bothDataSetsSummary %<>% mutate(metaP.neg.adj = p.adjust(metaP.neg, method="fdr"))
bothDataSetsSummary %>% arrange(metaP.neg.adj)
bothDataSets %>% filter(geneSymbol == "Hsp90ab1")
bothDataSetsSummary %<>% mutate(isSig = metaP.neg.adj < 0.05 | metaP.pos.adj < 0.05)
sum(bothDataSetsSummary$isSig)
nrow(bothDataSetsSummary)
bothDataSetsSummary %<>% mutate(direction = sign(metaP.neg.adj - metaP.pos.adj))

#add in medians across datasets
normMeans <- bothDataSets %>% select(geneSymbol, batch, Normal) %>% dcast(geneSymbol ~ batch ) %>% 
  rename(ChenBatch1FedAvg = B1, ChenBatch2FedAvg = B2, CampbellBatch6FedAvg = b6)
hungerMeans <- bothDataSets %>% select(geneSymbol, batch, Hungry) %>% dcast(geneSymbol ~ batch ) %>% 
  rename(ChenBatch1HungerAvg = B1, ChenBatch2HungerAvg = B2, CampbellBatch6HungerAvg = b6)
meanMatrix <- as_tibble(inner_join(normMeans, hungerMeans))
meanMatrix %<>% select(geneSymbol, noquote(order(colnames(meanMatrix))))
bothDataSetsSummary <- inner_join(meanMatrix, bothDataSetsSummary)

directionGroupingRatioAll <- bothDataSetsSummary %>% filter(isSig) %>% group_by(direction) %>% summarize(n=n()) %>% print()
directionGroupingRatioAll <- binom.test(directionGroupingRatioAll$n)$estimate

isSigGroupingRatioAll <- bothDataSetsSummary %>% group_by(isSig) %>% summarize(n=n()) %>% print()
isSigGroupingRatioAll <- binom.test(isSigGroupingRatioAll$n)$estimate

#fullSummaryTable <- 
#write out for genesets of interest
allGenes <- c()
for (targetGeneList in targetGeneLists) {
  
  targetSymbols <- read_csv(paste0("./data/genelists/",targetGeneList), col_names=F)$X1
  print(targetGeneList)  
  allGenes <- union(allGenes, targetSymbols)
  filteredRaw <- bothDataSets %>% filter(geneSymbol %in% targetSymbols)
  write_csv(filteredRaw, paste0("./results/cell types/", targetGeneList, ".hungerPerGene.wilcox.raw.combined.csv"))
  
  filteredSummary <- bothDataSetsSummary %>% filter(geneSymbol %in% targetSymbols)
  #direction bias test:
  directionGrouping <- filteredSummary %>% filter(isSig) %>% group_by(direction) %>% summarize(n=n()) 
  if (nrow(directionGrouping) ==2 ){
    print(paste("binomial test on direction p=", binom.test(directionGrouping$n, p=directionGroupingRatioAll )$p.value)) #value from all gene result
  }
  write_csv(filteredSummary, paste0("./results/cell types/", targetGeneList, ".hungerPerGene.wilcox.metaP.summary.csv"))
  
  isSigGrouping <- filteredSummary %>% group_by(isSig) %>% summarize(n=n()) 
  print(paste("Number genes isSig:", signif(nrow(filteredSummary %>% filter(isSig)))))
  print(paste("Percent genes isSig:", signif(nrow(filteredSummary %>% filter(isSig))/nrow(filteredSummary)*100, digits=4)))
  print(paste("binomial test on isSig p=", binom.test(isSigGrouping$n, p=isSigGroupingRatioAll )$p.value)) #value from all gene result
}

forSupplement <- bothDataSetsSummary %>% filter(geneSymbol %in% allGenes)
for (targetGeneList in targetGeneLists) {
  targetSymbols <- read_csv(paste0("./data/genelists/",targetGeneList), col_names=F)$X1
  targetGeneList <- gsub(".hypenFixed.mouse.txt", "", targetGeneList) 
  forSupplement %<>% mutate(!!targetGeneList := geneSymbol %in% targetSymbols)
}
write_csv(forSupplement, paste0("./results/cell types/SummarySupplement.csv"))

#############################
#specific code for Negraes - direcion testing
#deal with gene conversion - watch out for one to many
humanSymbolMap <- mouse2human(bothDataSets$geneSymbol) %>% rename(geneSymbol = mouseGene) %>% distinct()

targetGeneList <- "Negraes et al. Table S5.hypenFixed.mouse.txt"
targetSymbols <- read_csv(paste0("./data/genelists/",targetGeneList), col_names=F)$X1
filteredSummary <- bothDataSetsSummary %>% filter(geneSymbol %in% targetSymbols)

filteredSummary <- as_tibble(right_join(humanSymbolMap, filteredSummary))

#for gene expression studies, test direction against previous results
NegraesTable <- read_csv("./data/genelists/Negraes et al. Table S5.csv")
NegraesTable %<>% mutate(directionNegraes = sign(`Log2 Fold-Change`))
NegraesTable %<>% select(humanGene=Symbol, directionNegraes)
NegraesTable %<>% mutate(humanGene = gsub("-AS1$","", humanGene))
NegraesTable %<>% mutate(humanGene = gsub("-AS$","", humanGene))
NegraesTable %<>% mutate(humanGene = gsub("-AS2$","", humanGene))
NegraesTable %<>% mutate(humanGene = gsub("-IT1$","", humanGene))
filteredSummary <- left_join(filteredSummary, NegraesTable)

filteredSummary %<>% mutate(Negraes_validation = directionNegraes == direction)

filteredSummary %>% filter(isSig) %>% group_by(directionNegraes, direction) %>% summarize(n=n())

filteredSummary %>% filter(Negraes_validation, isSig) %>% select(geneSymbol, isSig, direction, directionNegraes)
trueHits <- nrow(filteredSummary %>% filter(Negraes_validation, isSig))
#manually ran a chisq.test test, not significant for agreement
