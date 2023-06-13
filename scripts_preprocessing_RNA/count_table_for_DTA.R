rm(list = ls())
library(tidyverse)

# load count table and change the colnames so they match to those in metadata
countTable = read.csv("countTables/countCerevisiae_genes_minusStranded.csv") %>% 
  column_to_rownames(var = "X")

changeColnames = unlist(strsplit(colnames(countTable), split = "_"))
changeColnames = changeColnames[seq(1, length(changeColnames), by = 2)]

colnames(countTable) = changeColnames

# load metadata and number of mapped rads
mappedReads = read.csv("mapping_stats/mappingStats_fullReport.csv")

metadata = read.csv("metadata.csv") %>% 
  full_join(mappedReads) %>% 
  mutate(normFactor = mappedReads_pombe/2000000)

rownames(metadata) = metadata$sampleNames

rm(list = setdiff(ls(), c("countTable", "metadata")))

countTable = countTable[, rownames(metadata)]

# normalize counts: 
normCount_full = sweep(countTable, 2, metadata$normFactor, "/")
write.csv(normCount_full, "countTables/countCerevisiae_genes_minusStranded_pombeNormalized.csv")

# let's have a look at the counts
row_means = rowMeans(countTable)
row_med = apply(countTable, 1, median)
numOfZeros = rowSums(countTable == 0)

plotCountOverview = data.frame(gene = rownames(countTable), med = row_med, means = row_means, zeros = numOfZeros)

p = ggplot(plotCountOverview, aes(med))
p + geom_histogram(bins = 100) +
  theme_bw() +
  xlim(-1, 10000)

p = ggplot(plotCountOverview, aes(means))
p + geom_histogram(bins = 100) +
  theme_bw() +
  xlim(-1, 25)

p = ggplot(plotCountOverview, aes(zeros))
p + geom_histogram(bins = 100) +
  theme_bw() +
  xlim(-1, 30)

# =========== filter out counts where we see 0 in more than half of the samples
countTableFilt = countTable[rowSums(countTable == 0) < 12,]
write.csv(countTableFilt, "countTables/filtered_countTables/countCerevisiae_genes_minusStranded__lessThan12zeros.csv")

# explore normalized counts of the filtered count table
normCount = sweep(countTableFilt, 2, metadata$normFactor, "/")
write.csv(normCount, "countTables/filtered_countTables/countCerevisiae_genes_minusStranded__lessThan12zeros_pombeNormalized.csv")

plot_normCounts = normCount %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(!gene, names_to = "sampleNames", values_to = "normCount") %>% 
  full_join(metadata)

p = ggplot(plot_normCounts, aes(x = sampleNames, y = normCount, color = sampleType))
p + geom_boxplot() +
  theme_bw() +
  ylim(0, 250) +
  facet_wrap(~thiouracil)
ggsave("figures/normCounts_boxplot_cropped250.png", width = 10, height = 3)


# plot replicates against each other
plotReplicates = normCount %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(!gene, names_to = "sampleNames", values_to = "normCount") %>%
  left_join(metadata %>% select(sampleNames, sampleType, timeInGal, thiouracil, replicates)) %>% 
  select(-sampleNames) %>% 
  pivot_wider(names_from = replicates, values_from = normCount) %>% 
  mutate(lab = paste0(sampleType, "_", timeInGal, "min"))
p = ggplot(plotReplicates, aes(x = rep1, y = rep2, color = sampleType)) 
p + geom_point()+
  theme_bw() +
  facet_grid(thiouracil~lab) +
  xlim(0, 10000) +
  ylim(0, 10000)
ggsave("figures/normCounts_compareReplicates.png", width = 10, height = 4)


##### create count table only for samples with thiouracil added 
rm(list = setdiff(ls(), c("countTable", "metadata")))

metadataThio = metadata %>% 
  filter(thiouracil == "yes")
filtThio = countTable[, rownames(metadataThio)]

# explore normalized counts of the filtered count table
normCount = sweep(filtThio, 2, metadataThio$normFactor, "/")
write.csv(normCount, "countTables/filtered_countTables/countCerevisiae_genes_minusStranded__ThiouracilYESonly_pombeNormalized.csv")



