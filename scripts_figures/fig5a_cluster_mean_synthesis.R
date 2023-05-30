rm(list = ls())

library(tidyverse)
library(viridis)
# load table with assignement of genes to heatmap clusters: from script "fig5a_heatmap.R"
heatmapClusters = read.csv("data/analysis/heatmap_clusters/heatmapClusters.csv")
# get mean synthesis rates for each cluster
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin)) %>% 
  inner_join(heatmapClusters) %>% 
  group_by(cluster) %>% 
  summarise(meanSynthesis = mean(cropped)) %>% 
  mutate(clusterName = paste0("cluster ", cluster))

p = ggplot(synthesisRates,aes(y = reorder(cluster, -cluster), x = meanSynthesis, fill = meanSynthesis))
p + geom_bar(stat = "identity") +
  scale_fill_distiller(palette = "Reds", direction = 1)+
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9)) +
  ylab("cluster number") +
  xlab(paste("mean synthesis rate", "[mRNA/cell/min]", sep = "\n"))
#ggsave("figures/panels/fig6/meanSynthesis.pdf", width = 6, height = 8, units = "cm")
