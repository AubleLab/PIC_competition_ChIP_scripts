rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)


set.seed(42)


# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-peakName) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum)


###################################################
#------ heatmap ------
###################################################
heatmapTable = plotResTime %>%
  ungroup() %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  column_to_rownames(var = "gene")
naPerRow = rowSums(is.na(heatmapTable))

# filter rows with more than 1 NA
heatmapTableFilt = heatmapTable[naPerRow <= 1,c("TBP", "TFIIA", "TFIIB", "TFIIF", "TFIIE")]
genes = rownames(heatmapTableFilt)
rownames(heatmapTableFilt) = NULL

# plot - ward clustering showst the best results
Heatmap(heatmapTableFilt, cluster_columns = F, clustering_method_rows = "ward.D", name = "resTime")
# pdf(file = "/Users/lilcrusher/competitionChIP/figures/heatmap_all_GTF/heatmap_max1NA_wardClust.pdf", width = 4, height = 6)
# draw(Heatmap(heatmapTableFilt, cluster_columns = F, clustering_method_rows = "ward.D", name = "resTime"))
# dev.off()


# plot - ward clustering shows the best results + k-means clustering
Heatmap(heatmapTableFilt, 
        cluster_columns = F, 
        clustering_method_rows = "ward.D", 
        name = "resTime", 
        row_split = 9,
        row_gap = unit(3, "mm"), border = T)
# pdf(file = "/Users/lilcrusher/competitionChIP/figures/heatmap_all_GTF/heatmap_max1NA_wardClust_5clusters.pdf", width = 4, height = 6)
# draw(Heatmap(heatmapTableFilt, 
#              cluster_columns = F, 
#              clustering_method_rows = "ward.D", 
#              name = "resTime", 
#              row_split = 10,
#              row_gap = unit(3, "mm"), border = T))
# dev.off()

# extract the gene names for each cluster
ht = Heatmap(heatmapTableFilt, 
             cluster_columns = F, 
             clustering_method_rows = "ward.D", 
             name = "resTime", 
             row_split = 10,
             row_gap = unit(3, "mm"), border = T)
ht = draw(ht)
clusters = row_order(ht)
# # extract the gene names falling into each cluster
# for (i in 1:length(clusters)){
#   clusterRows = clusters[[i]]
#   clustGenes = genes[clusterRows]
#   
#   write.table(clustGenes, paste0("data/heatmap_clusters/cluster_", i, ".txt"), sep = "\n", 
#               row.names = F, col.names = F, quote = F)
# }
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
  write.csv(clustTableFull, "data/heatmapClusters.csv", row.names = F, quote = F)
}

plotTable = heatmapTableFilt %>% 
  mutate(gene = genes) %>% 
  pivot_longer(!gene, names_to = "factorName", values_to = "resTime", values_drop_na = T) %>% 
  left_join(clustTableFull) 

p = ggplot(plotTable, aes(x = factorName, y = resTime, fill = factorName))
p + geom_boxplot(outlier.shape=NA, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey70") + 
  theme_classic() +
  xlab(" ") +
  ylab("residence time [min]") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_wrap(~cluster, scales = "free", ncol = 1) +
  ylim(0,20)





