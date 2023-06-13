rm(list = ls())
library(tidyverse)

# load count table and change the colnames so they match to those in metadata
countTable = read.csv("data/RNA_count_tables/countCerevisiae_genes_minusStranded.csv") %>% 
  column_to_rownames(var = "X")

changeColnames = unlist(strsplit(colnames(countTable), split = "_"))
changeColnames = changeColnames[seq(1, length(changeColnames), by = 2)]

colnames(countTable) = changeColnames

# load metadata
metadata = read.csv("data/metadata.csv") %>% 
  mutate(normFactor = mappedReads_pombe/2000000)

rownames(metadata) = metadata$sampleNames

rm(list = setdiff(ls(), c("countTable", "metadata")))

# keep only samples which had thiouracil added
metadataThio = metadata %>% 
  dplyr::filter(added_4sU == "+")
filtThio = countTable[, rownames(metadataThio)]

# explore normalized counts of the filtered count table
normCount = sweep(filtThio, 2, metadataThio$normFactor, "/")
write.csv(normCount, "data/RNA_count_tables/countCerevisiae_genes_minusStranded__ThiouracilYESonly_pombeNormalized.csv")



