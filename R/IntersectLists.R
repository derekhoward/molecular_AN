library(gplots)

LutterBinge <- read_csv("./data/genelists/Lutter et al. Table S4.BingeEating.hypenFixed.txt", col_names=F)$X1
LutterRestricted <- read_csv("./data/genelists/Lutter et al. Table S3.RestrictedEating.hypenFixed.txt", col_names=F)$X1
Duncan <- read_csv("./data/genelists/Duncan et al. rs4622308.hypenFixed.txt", col_names=F)$X1
Negraes <- read_csv("./data/genelists/Negraes et al. Table S5.hypenFixed.txt", col_names=F)$X1
GSE60190 <- read_csv("./data/genelists/GSE60190.reanalysis.txt", col_names=F)$X1

venn( list(LutterRestricted=LutterRestricted, Duncan=Duncan, Negraes=Negraes, GSE60190 = GSE60190) )

combined <- union(union(union(LutterRestricted, Duncan),Negraes),GSE60190)

write_csv(tbl_df(sort(unique(combined))), "./data/genelists/Four.combined.txt", col_names = F)



