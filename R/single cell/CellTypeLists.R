library(dplyr)
library(readr)
library(homologene)
library(tmod)


#targetGeneList <- "Negraes et al. Table S5.hypenFixed.txt"
#targetGeneList <- "Duncan et al. rs4622308.hypenFixed.txt"
#targetGeneList <- "Lutter et al. Table S3.RestrictedEating.hypenFixed.txt"
#targetGeneList <- "Lutter et al. Table S4.BingeEating.hypenFixed.txt"

targetSymbolsHuman <- read_csv(paste0("./data/genelists/",targetGeneList), col_names=F)$X1


tmodNames <- data.frame()
modules2genes <- list()
otherGeneListsFolder <- "./data/single_cell/gene lists/"

#background set is human genes in homologene
length(background <- unique(tbl_df(homologeneData) %>% filter(Taxonomy == 9606) %>% .$Gene.Symbol))

for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "PH", geneListFilename) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    print(" converting from mouse to human")
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }
  
  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

result <- tmodHGtest(fg = targetSymbolsHuman, bg = background, mset = geneSets, qval = 1.01, filter = F)
(result <- tbl_df(result) %>% dplyr::select(Title, geneCount =B, overlap = b, P.Value, adj.P.Val))
result

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)

writeTableWrapper <- function(prefixFilter, result) {
  subsetResult <- filter(result, grepl(prefixFilter, Title))
  subsetResult$oldTitle <- subsetResult$Title
  subsetResult$Title <- gsub(paste0(prefixFilter, "."), "", subsetResult$Title)
  subsetResult$Title <- gsub("[.]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("[_]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("Neuron interneuron", "Interneuron", subsetResult$Title)
  subsetResult$Title <- gsub("ligo", "ligodendrocyte", subsetResult$Title)
  subsetResult$adj.P.Val <- p.adjust(subsetResult$P.Value, method="fdr")
  subsetResult$adj.P.Val <- signif(subsetResult$adj.P.Val, digits=3)
  #write_csv(dplyr::select(subsetResult, `Cell-type or class` = Title,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Val), paste0(baseFilename,".", prefixFilter,".csv")) 
  subsetResult
}

(writeTableWrapper("NeuroExpresso.All", result))
(writeTableWrapper("NeuroExpresso.Cortex", result))
(writeTableWrapper("NeuroExpresso.Amygdala", result))
(writeTableWrapper("NeuroExpresso.Brainstem", result))

sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Brainstem.Microglia, targetSymbolsHuman)) 
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Brainstem.Microglia_activation, targetSymbolsHuman)) #For Lutter
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Brainstem.Microglia_deactivation, targetSymbolsHuman)) #For Lutter
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Cortex.OligoPrecursors, targetSymbolsHuman)) #For Duncan
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Brainstem.Microglia, targetSymbolsHuman))
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Midbrain.Microglia, targetSymbolsHuman))
sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Brainstem.Noradrenergic, targetSymbolsHuman))

sort(intersect(geneSets$MODULES2GENES$NeuroExpresso.Cortex.Oligo, targetSymbolsHuman))

