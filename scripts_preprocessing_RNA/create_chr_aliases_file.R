rm(list = ls())

library(tidyverse)


GTF = read.delim("/Users/lilcrusher/annotations/yeast/Saccharomyces_cerevisiae.R64-1-1.106.gtf",skip = 5, header = F)

chromosomes = GTF %>% 
  select(V1) %>% 
  distinct() %>% 
  arrange(V1)


chromosomeManual = c("I", "II", "III", "IV", "V", 
                     "VI", "VII", "VIII", "IX", "X", 
                     "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito")


chrTable = data.frame(chromosomeManual, chromosomeManual)

write.table(chrTable, "/Users/lilcrusher/annotations/yeast/chrAliasesForFeatureCounts.txt", sep = ",", quote = F, 
            row.names = F,col.names = F)
