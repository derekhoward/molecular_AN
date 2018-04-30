library(dplyr)
library(readr)
library(homologene)

genelists <- list.files("./data/genelists/", pattern= ".*genelist.*", full.names = T)

for(genelist in genelists) {
  
  baseName <- gsub("[.]txt", "", genelist)
  baseName <- gsub("[.]genelist", "", baseName)
  print(baseName)
  #replace -AS, -AS1 -AS2 and -IT1 with nothing
  targetSymbolsHuman <- read_csv(genelist, col_names=F)$X1
  targetSymbolsHuman <- gsub("-AS1$","", targetSymbolsHuman)
  targetSymbolsHuman <- gsub("-AS$","", targetSymbolsHuman)
  targetSymbolsHuman <- gsub("-AS2$","", targetSymbolsHuman)
  targetSymbolsHuman <- gsub("-IT1$","", targetSymbolsHuman)
  targetSymbolsHuman <- unique(targetSymbolsHuman)
  targetSymbolsHumanOrig <- targetSymbolsHuman
  write.table(sort(unique(targetSymbolsHumanOrig)), paste0(baseName, ".hypenFixed.txt"), col.names = F, row.names = F, quote = F)
  
  orthologs <- tbl_df(human2mouse(targetSymbolsHuman))
  setdiff(targetSymbolsHuman,orthologs$humanGene)
  targetSymbolsMouse <- unique(orthologs$mouseGene)
  write.table(sort(unique(targetSymbolsMouse)), paste0(baseName, ".hypenFixed.mouse.txt"), col.names = F, row.names = F, quote = F)
}

