rm(list = ls())

library(tidyverse)

Lickwar = read.csv("data/Rap1_resTime.csv") %>% 
  mutate(gene = ifelse(startsWith(gene,"i"), sub(".", "", gene), gene)) %>% 
  dplyr::select(gene, resTime) %>% 
  dplyr::rename(Lickwar = resTime)

Rap1 = read.csv("data/Fkh1_Rap1_resTimes.csv") %>% 
  dplyr::filter(factorName == "Rap1")


plotTable = Rap1 %>% 
  inner_join(Lickwar)


r = cor(plotTable$resTime, plotTable$Lickwar)
rho = cor(plotTable$resTime, plotTable$Lickwar, method = "spearman")


p = ggplot(plotTable, aes(x = resTime, y = Lickwar))
p + geom_point() +
  theme_classic() +
  ggtitle(paste0("r = ", round(r, digits = 2)))
ggsave("figures/panels/random/compareRap1.pdf", width = 3, height = 3)
