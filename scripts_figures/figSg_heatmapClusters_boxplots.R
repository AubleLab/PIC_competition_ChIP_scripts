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


###################################################
#------ heatmap ------
###################################################
heatmapTable = plotResTime %>%
  ungroup() %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  column_to_rownames(var = "gene")
naPerRow = rowSums(is.na(heatmapTable))

# filter rows with NAs in them
heatmapTableFilt = heatmapTable[naPerRow == 0,c("TFIIA", "TFIIB", "TFIIF", "TFIIE")]
genes = rownames(heatmapTableFilt)
rownames(heatmapTableFilt) = NULL

# z-score normalize the values
zHeatmap = scale(heatmapTableFilt)

col_fun = colorRamp2(c(-2, 0, 2), c("#31A9B8", "white", "#CF3721"))
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
}

plotTable = as.data.frame(zHeatmap) %>% 
  mutate(gene = genes) %>% 
  pivot_longer(!gene, names_to = "factorName", values_to = "resTime", values_drop_na = T) %>% 
  left_join(clustTableFull) %>% 
  mutate(cluster = paste0("cluster ", cluster))
plotTable$factorName = factor(plotTable$factorName, levels = c("TFIIA", "TFIIB", "TFIIF", "TFIIE"))
plotTable$cluster = factor(plotTable$cluster, levels = paste0(rep("cluster ", 10), seq(1, 10)))

p = ggplot(plotTable, aes(x = factorName, y = resTime, color = factorName, fill = factorName))
p + geom_boxplot(outlier.shape=NA, alpha = 0.2) +
  theme_classic() +
  xlab(" ") +
  ylab("z(residence time)") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_wrap(~cluster, scales = "free", ncol = 5) + 
  scale_color_manual(values = c("#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")) + 
  scale_fill_manual(values = c("#E6772E", "#3D4C53", "#4DB3B3", "#E64A45"))
ggsave("figures/panels/figSf/heatmapClusters_boxplots.pdf", width = 17, height = 8, units = )




