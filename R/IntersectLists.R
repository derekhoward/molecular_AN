library(gplots)

LutterBinge <- read_csv("./data/genelists/Lutter et al. Table S4.BingeEating.hypenFixed.txt", col_names=F)$X1
LutterRestricted <- read_csv("./data/genelists/Lutter et al. Table S3.RestrictedEating.hypenFixed.txt", col_names=F)$X1
Duncan <- read_csv("./data/genelists/Duncan et al. rs4622308.hypenFixed.txt", col_names=F)$X1
Negraes <- read_csv("./data/genelists/Negraes et al. Table S5.hypenFixed.txt", col_names=F)$X1
Watson <- read_csv("./data/genelists/Watson et al.TableS6.protein_coding.ALL.hypenFixed.txt", col_names=F)$X1

intersect(Negraes, LutterRestricted)
intersect(Negraes, LutterBinge)
intersect(Negraes, Watson)

venn( list(LutterRestricted=LutterRestricted, Duncan=Duncan, Negraes=Negraes, Watson=Watson) )
venn( list(LutterRestricted=LutterRestricted, Duncan=Duncan, Negraes=Negraes, LutterBinge=LutterBinge) )

combined <- union(union(union(LutterRestricted, Duncan),Negraes),Watson)

write_csv(tbl_df(sort(unique(combined))), "./data/genelists/Four.combined.txt", col_names = F)


length(Negraes)
length(LutterBinge)
length(LutterRestricted)
length(Watson)
