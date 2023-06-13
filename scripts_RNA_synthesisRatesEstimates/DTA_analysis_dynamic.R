rm(list = ls())

library(DTA)
library(tidyverse)

# upload Pombe normalized count table with samples, where thiouracil was added
countTable = read.csv("data/RNA_count_tables/countCerevisiae_genes_minusStranded__ThiouracilYESonly_pombeNormalized.csv") %>% 
  column_to_rownames(var = "X")

# load metadata and filter out the samples where thiouracil was not added
metadata = read.csv("data/metadata.csv") %>% 
  filter(added_4sU == "+") %>% 
  dplyr::rename(name = sampleNames) %>% 
  mutate(time = 6) %>% 
  mutate(timeframe = paste0("00-", as.character(timeInGal))) %>% 
  mutate(timecourse = ifelse(timeInGal == 20, 1, 2))
rownames(metadata) = metadata$name
metadata$nr = rep(c(1,2,3,4), 3)

# there cannot be any zero counts in the matrix - let's filter those genes out
zeroSum = rowSums(countTable == 0)
zeroSum = zeroSum[zeroSum == 0]

countTable_nonZero = countTable[names(zeroSum),]

# create list of reliable genes - all genes left
# reformat metadata and count table to fit the DTA function
metadata_DTA = as.matrix(metadata)
reliable_DTA_plusOne = rownames(countTable_nonZero)
countTable_DTA_plusOne = as.matrix(countTable_nonZero)

res_DTA = DTA.dynamic.estimate(phenomat = metadata_DTA,
                       datamat = countTable_DTA_plusOne,
                       tnumber = Sc.tnumber, 
                       check = T,
                       ccl = 150,mRNAs = 60000,
                       reliable = reliable_DTA_plusOne,
                       condition = "real_data",save.plots = TRUE,
                       notinR = TRUE,folder = "data/figures", 
                       ratiomethod = "bias")
# extract 20 min results
resList = res_DTA$`1`

# extract gene names and synthesis rates (divide by 150 to make the result
# per minute not per cell cycle)
synthesis = resList$sr
genes = names(synthesis)

resTable20 = data.frame(gene = genes,
                        synthesis20 = synthesis/150)

# extract 60 min results
resList = res_DTA$`2`

synthesis = resList$sr
genes = names(synthesis)

resTable = data.frame(gene = genes, 
                      synthesis60 = synthesis/150) %>% 
  full_join(resTable20)

# write.csv(resTable, "data/synthesisRates_60_vs_20_min.csv")
