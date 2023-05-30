rm(list = ls())
library(tidyverse)

# upload gene classes

geneClass = read.csv("data/analysis/heatmap_clusters/TF_enrichment/sigP_TFenrichment_labels.csv")
# plot results from my TF enrichment
sigRes = read.csv("data/analysis/heatmap_clusters/TF_enrichment/sigP_TFenrichment.csv")
plotRes = sigRes %>% 
  dplyr::select(TF, p, cluster, logFDR, FDR) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(TF = str_to_title(TF)) %>% 
  mutate(logP = -log10(p)) %>% 
  mutate(sigFDR = ifelse(FDR < 0.05, "sig", "nonSig")) %>% 
  left_join(geneClass) %>% 
  drop_na() %>% 
  mutate(combo = paste(cluster, label1, TF))
  

p = ggplot(plotRes, 
           aes(x = combo, 
               y = 1, fill = logFDR))
p + geom_bar(stat= "identity", color = "grey80") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_distiller(palette = "GnBu", direction = 1, na.value = "grey80",
                       name = expression(paste(-log[10],"(FDR)"))) +
  ylab("")+
  xlab("") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        axis.text=element_text(size=8), legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(cols = vars(cluster), rows = vars(label1), scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("grey50", "white")) +
  scale_x_discrete(breaks=plotRes$combo,
                     labels=plotRes$TF) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 1))
ggsave("figures/panels/figSf/myTFenrichment_clusters.pdf", width = 16, height = 10, units = "cm")


