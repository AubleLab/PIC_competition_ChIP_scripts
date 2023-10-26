rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)


set.seed(42)


# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-c(peakName, TBP)) %>% 
  pivot_longer(cols = TFIIA:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum)

# upload synthesis rates
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin)) %>% 
  right_join(plotResTime) %>% 
  dplyr::select(gene, cropped) %>% 
  distinct()
rownames(synthesisRates) = synthesisRates$gene
###################################################
#------ heatmap ------
###################################################
heatmapTable = plotResTime %>%
  ungroup() %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  column_to_rownames(var = "gene")
naPerRow = rowSums(is.na(heatmapTable))
# make sure that synthesis rate is the same order
synthesisOrdered = synthesisRates[rownames(heatmapTable),]

# filter rows with NAs in them
heatmapTableFilt = heatmapTable[naPerRow == 0,c("TFIIA", "TFIIB", "TFIIF", "TFIIE")]
genes = rownames(heatmapTableFilt)
rownames(heatmapTableFilt) = NULL

synthesisPlot = synthesisRates[naPerRow == 0,"cropped"]

# z-score normalize the values
zHeatmap = scale(heatmapTableFilt)

col_fun = colorRamp2(c(-2, 0, 2), c("#31A9B8", "white", "#CF3721"))

synthesis_col = colorRamp2(c(min(synthesisPlot, na.rm = TRUE),
                             (min(synthesisPlot, na.rm = TRUE) + (max(synthesisPlot, na.rm = TRUE)))/2, 
                             max(synthesisPlot, na.rm = TRUE)), 
                     magma(n = 5, direction = -1)[c(1,3,5)])
row_ha = rowAnnotation(synthesis = synthesisPlot, col = list(synthesis = synthesis_col))
# plot - ward clustering shows the best results + k-means clustering

Heatmap(zHeatmap, 
        cluster_columns = F, 
        clustering_method_rows = "ward.D", 
        name = "z(res. time)", 
        row_split = 10,
        row_gap = unit(1.5, "mm"), border = T, col = col_fun)#,
        #right_annotation = row_ha)


# extract the gene names for each cluster
ht = Heatmap(zHeatmap, 
             cluster_columns = F, 
             clustering_method_rows = "ward.D", 
             name = "z(res. time)", 
             row_split = 10,
             row_gap = unit(1, "mm"), border = T, col = col_fun)
ht = draw(ht)
clusters = row_order(ht)
# extract the gene names falling into each cluster
for (i in 1:length(clusters)){
  clusterRows = clusters[[i]]
  clustGenes = genes[clusterRows]
  clustTable = data.frame(cluster = rep(i, length(clustGenes)) ,gene = clustGenes)
  if (i == 1){
    clustTableFull = clustTable
  } else {
    clustTableFull = rbind(clustTableFull, clustTable)
  }
  write.table(clustGenes, paste0("data/analysis/heatmap_clusters/", i, ".txt"), quote = F, row.names = F, col.names = F)
}
write.csv(clustTableFull, "data/analysis/heatmap_clusters/heatmapClusters.csv", row.names = F, quote = F)


