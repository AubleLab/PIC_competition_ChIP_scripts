rm(list = ls())

library(tidyverse)


resTable = read.csv("data/residence_times_all.csv") %>% 
  select(-peakName) %>% 
  pivot_longer(TBP:TFIIF, names_to = "TF", values_to = "resTime") %>% 
  drop_na(resTime) %>% 
  distinct() %>% 
  pivot_wider(names_from = TF, values_from = resTime) %>% 
  replace(is.na(.), " ")
write.csv(resTable,"manuscript/Supplements/Table_S1.csv", quote = F, row.names = F)
