rm(list = ls())
library(tidyverse)

# plot results from my TF enrichment
# upload annotation (manually created labels)
annot = read.csv("data/analysis/heatmap_clusters/TF_enrichment/sigP_TFenrichment_labels.csv")
# upload results of significant TF target enrichment (Yeats Epigeneome database gene targets): 
#these were generated with "fig5_heatmapClusters_TFenrichment.R"
sigRes = read.csv("data/analysis/heatmap_clusters/TF_enrichment/sigP_TFenrichment.csv")
plotRes = sigRes %>% 
  dplyr::select(TF, p, cluster, logFDR, FDR) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(TF = str_to_title(TF)) %>% 
  mutate(logP = -log10(p)) %>% 
  mutate(sigFDR = ifelse(FDR < 0.05, "sig", "nonSig")) %>% 
  inner_join(annot) %>% 
  dplyr::filter(label1 != "Pol II" & !startsWith(label1, "TF") & label1 != "TBP")

p = ggplot(plotRes, 
           aes(y = reorder(TF, logP), 
                        x = logP, fill = logFDR, color = sigFDR))
p + geom_bar(stat= "identity") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_distiller(palette = "GnBu", direction = 1, na.value = "grey80",
                       name = expression(paste(-log[10],"(FDR)"))) +
  ylab("")+
  xlab(expression(paste(-log[10],"(p)"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        axis.text=element_text(size=8), legend.position = "top") + 
  facet_grid(rows = vars(cluster), scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("grey50", "white"))
#ggsave("figures/panels/fig5/myTFenrichment_clusters.pdf", width = 5, height = 10, units = "cm")
#ggsave("figures/panels/fig5/myTFenrichment_clusters_legend.pdf", width = 10, height = 10, units = "cm")

