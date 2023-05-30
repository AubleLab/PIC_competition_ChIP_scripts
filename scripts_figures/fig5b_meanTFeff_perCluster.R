rm(list = ls())

library(tidyverse)
library(viridis)

set.seed(42)

# load synthesis rate
synthesisRates = read.csv("data/synthesisRates.csv") 

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-TBP) %>% 
  pivot_longer(cols = TFIIA:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum)

# upload heatmap clusters
heatmapClusters = read.csv("data/analysis/heatmap_clusters/heatmapClusters.csv")

# upload synthesis rates, attach residence times, 
# compute TR efficiency as synthesis rate*residence time, 
# get mean and median TR efficiency 
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  inner_join(heatmapClusters) %>% 
  inner_join(plotResTime) %>% 
  mutate(TR_eff = resTimesNum * synthesis_perCell_perMin) %>% 
  group_by(cluster, factorName) %>% 
  summarise(meanEfficiency = mean(TR_eff),
            medEfficiency = median(TR_eff)) %>% 
  mutate(clusterName = paste0("cluster ", cluster))

synthesisRates$factorName = factor(synthesisRates$factorName, levels = c("TFIIA", "TFIIB", "TFIIF", "TFIIE"))

# plot
p = ggplot(synthesisRates,aes(y = reorder(cluster, -cluster), x = medEfficiency, fill = meanEfficiency))
p + geom_bar(stat = "identity", color = "grey80") +
  scale_fill_distiller(palette = "Blues", direction = 1)+
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        legend.position = "none") +
  ylab("cluster number") +
  xlab(paste("median transcription efficiency", "[mRNA per binding event]", sep = "\n")) +
  labs(fill = paste("mean transcription efficiency", "[mRNA per binding event]", sep = "\n")) +
  facet_wrap(~factorName)
#ggsave("figures/panels/fig5/meanTRefficiency_perCluster.pdf", width = 7, height = 8, units = "cm")

