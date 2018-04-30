#requires manually writing out single cell data for the gene
library(ggsignif)
library(ggplot2)

Rps26Chen <- read_csv("./results/cell types/Rps26.GSE87544.csv")

Rps26Chen %>% group_by(batch, group_combined) %>% summarize(n=n(), mean = mean(log1Expression), median = median(log1Expression)) 
Rps26Chen %>% group_by(batch, group_combined) %>% summarize(n=n())

#compare batch 3 and 4
Rps26Chen %>% ungroup() %>% filter(batch %in% c("B3", "B4")) %>% summarize(
  p = wilcox.test(log1Expression ~ batch)$p.value,
)
Rps26Chen %>% filter(batch %in% c("B3", "B4")) %>% group_by(batch) %>% summarize(
  mean = mean(log1Expression),
  group_combined = first(group_combined)
)

Rps26ChenSummary <- Rps26Chen %>% group_by(batch, group_combined) %>% filter(!batch %in% c("b1", "b2", "b3")) %>% summarize(
  p = wilcox.test(log1Expression ~ group_combined)$p.value,
)

wilcox.test(log1Expression ~ group_combined, Rps26Chen)
t.test(log1Expression ~ group_combined, Rps26Chen)

#Campbell
Rps26Campbell <- read_csv("./results/cell types/Rps26.GSE93374.csv")
Rps26Campbell %>% filter(batches=="b6", sex=="F")
hist(Rps26Campbell$log1Expression)

Rps26Campbell %>% group_by(batches, Sex_pred, FvF) %>% summarize( n= n())
Rps26Campbell %>% group_by(batches, sex, FvF) %>% summarize( n= n())
Rps26Campbell %>% group_by(batches, sex, group) %>% summarize( n= n())

boxplot(log1Expression ~ FvF,Rps26Campbell)

Rps26CampbellSummary <- Rps26Campbell %>% group_by(batches, Sex_pred) %>% filter(!batches %in% c("b1", "b2", "b3")) %>% summarize(
  p = wilcox.test(log1Expression ~ FvF)$p.value,
) 

meansPerGroup <- Rps26Campbell %>% group_by(FvF, Sex_pred, batches) %>% summarize(
  meanExp = mean(log1Expression)
)
meansPerGroup %<>% spread( FvF, meanExp)
Rps26CampbellSummary <- inner_join(meansPerGroup, Rps26CampbellSummary) %>% print()

#make a violin plot
forPlot <- bind_rows(
  Rps26Campbell %>% filter(batches == "b6", sex == "F") %>% select(log1Expression, condition = FvF) %>% mutate(study = "Campbell, batch 6")
  ,
  Rps26Chen %>% filter(batch == "B2")  %>% select(log1Expression, condition = group_combined) %>% mutate(study = "Chen, batch 2")
  ,
  Rps26Chen %>% filter(batch == "B1")  %>% select(log1Expression, condition = group_combined) %>% mutate(study = "Chen, batch 1")
)

forPlot %<>% mutate(condition = if_else(condition %in% c("Fast", "Hungry"), "food deprived", condition))
forPlot %<>% mutate(condition = if_else(condition %in% c("Fed", "Normal"), "ad libitum", condition))
forPlot$study <- factor(forPlot$study, levels=rev(unique(forPlot$study)))

ggplot(forPlot, aes(y=log1Expression, x=condition, color=condition)) + 
  geom_violin(adjust = 1/2, draw_quantiles = c(0.5)) + 
  theme_bw() + ylab("Rps26 Expression (log)") + theme(legend.position="none") + xlab("") +
  geom_signif(comparisons = list(unique(forPlot$condition)), color="black", test='wilcox.test', textsize = 3.2) +
  facet_wrap(~study, scales = "free_y") 
ggsave("/Users/lfrench/Google Drive/manuscripts/MolecularAnorexia/figures/Rps26Focus.pdf", plot = last_plot(), device = "pdf", width = 6, height = 5)
#save as 6x5 PDF
