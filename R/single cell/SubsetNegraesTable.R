library(homologene)

#spits Negraes table into up and down regulated
colnames(NegraesTable <- read_csv("./data/genelists/Negraes et al. Table S5.csv"))
NegraesTable %<>% mutate(direction = sign(`Log2 Fold-Change`))

NegraesTable %<>% mutate(BaseSymbol = Symbol)

NegraesTable %<>% mutate(BaseSymbol =  gsub("-AS1$","", BaseSymbol))
NegraesTable %<>% mutate(BaseSymbol =  gsub("-AS$","", BaseSymbol))
NegraesTable %<>% mutate(BaseSymbol =  gsub("-AS2$","", BaseSymbol))
NegraesTable %<>% mutate(BaseSymbol =   gsub("-IT1$","", BaseSymbol))

NegraesTable %<>% select(`Ensembl ID`, Symbol, BaseSymbol, everything())

#mapping <- human2mouse(NegraesTable$BaseSymbol)
#mapping %<>% group_by(humanGene) %>% summarize(mouseGenes = list(mouseGene)) %>% rename(BaseSymbol = humanGene)
#NegraesTable <- left_join(NegraesTable, mapping)

#mark microglia genes
microgliaMarkers <- read.csv("./data/single_cell/gene lists/NeuroExpresso.Brainstem.Microglia.txt",header=F,stringsAsFactors = F)
microgliaMarkers <-mouse2human(microgliaMarkers$V1)$humanGene

microgliaMarkersActivated <- read.csv("./data/single_cell/gene lists/NeuroExpresso.Brainstem.Microglia_activation.txt",header=F,stringsAsFactors = F)
microgliaMarkersActivated <- mouse2human(microgliaMarkersActivated$V1)$humanGene

microgliaMarkersDeactivated <- read.csv("./data/single_cell/gene lists/NeuroExpresso.Brainstem.Microglia_deactivation.txt",header=F,stringsAsFactors = F)
microgliaMarkersDeactivated <- mouse2human(microgliaMarkersDeactivated$V1)$humanGene

NegraesTable %<>% mutate(isMicrogliaMarker = BaseSymbol %in% microgliaMarkers)
NegraesTable %<>% mutate(isActivatedMicrogliaMarker = BaseSymbol %in% microgliaMarkersActivated)
NegraesTable %<>% mutate(isDeactivatedMicrogliaMarker = BaseSymbol %in% microgliaMarkersDeactivated)
NegraesTable %<>% filter(isMicrogliaMarker | isActivatedMicrogliaMarker | isDeactivatedMicrogliaMarker) %>% arrange(desc(isMicrogliaMarker), isActivatedMicrogliaMarker)
write_csv(NegraesTable, "./results/Negraes et al. Table S5.microglia.csv") 

