rm(list = ls())
library(tidyverse)

# load results obtaiend from g:Profiler (for genes from each heatmap cluster)
enrichmentFiles = list.files("data/analysis/heatmap_clusters/gProfiler", "*.csv", full.names = TRUE)

for (i in enrichmentFiles){
  clusterNum = unlist(strsplit(basename(i), split = ".", fixed = TRUE))[1]
  
  
  if (i == enrichmentFiles[1]){
    resTable = read.csv(i) %>% 
      mutate(cluster = clusterNum)
  } else {
    
    helpTable = read.csv(i) %>% 
      mutate(cluster = clusterNum)
    
    resTable = rbind(resTable, helpTable)
  }

}

gprofiler =  resTable %>% 
  dplyr::filter(source == "TF") %>% 
  separate(term_name, into = c("toss", "TF"), sep = ": ", extra = "drop") %>% 
  dplyr::select(-toss) %>% 
  separate(TF, into = c("TF", "toss"), sep = ";", extra = "drop") %>% 
  dplyr::select(-toss) %>% 
  mutate(TF = ifelse(endsWith(TF, "p"), str_sub(TF, end = -2), TF)) %>% 
  group_by(TF, cluster) %>% 
  dplyr::slice(which.max(negative_log10_of_adjusted_p_value)) %>% 
  dplyr::select(TF, negative_log10_of_adjusted_p_value, cluster)%>% 
  mutate(TF = str_to_title(TF))



# compare results from our TF enrichment and from g:Profiler
# these were generated with "fig5_heatmapClusters_TFenrichment.R"
sigRes = read.csv("data/analysis/heatmap_clusters/TF_enrichment/sigTFenrichment.csv")
plotRes = sigRes %>% 
  dplyr::select(TF, logFDR, cluster) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(TF = str_to_title(TF)) %>% 
  full_join(gprofiler) %>% 
  replace_na(list(negative_log10_of_adjusted_p_value = 0, logFDR = 0))

# plot results from g:Profiler validated with our results
plotRes = sigRes %>% 
  dplyr::select(TF, logFDR, cluster) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(TF = str_to_title(TF)) %>% 
  mutate(TF = ifelse(TF == "Matalpha2" | TF == "Mcm1", "Matalpha2-Mcm1", TF)) %>% 
  full_join(gprofiler) %>% 
  dplyr::filter(!is.na(negative_log10_of_adjusted_p_value))
p = ggplot(plotRes, aes(y = reorder(TF, negative_log10_of_adjusted_p_value), 
                        x = negative_log10_of_adjusted_p_value, fill = logFDR))
p + geom_bar(stat= "identity", color = "grey80", linewidth = 1) +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_distiller(palette = "GnBu", direction = 1, na.value = "grey80",
                       name = paste("validation", "-log10(FDR)", sep = "\n")) +
  ylab("")+
  xlab(expression(paste(-log[10],"(padj)"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        axis.text=element_text(size=8), legend.position = "top") + 
  facet_grid(rows = vars(cluster), scales = "free_y", space = "free_y")

# upload KEGG and WP results from g:Profiler
pathways = resTable%>% 
  dplyr::filter(source %in% c("GO:BP", "GO:MF", "KEGG", "WP")) %>% 
  group_by(source, cluster) %>% 
  arrange(desc(negative_log10_of_adjusted_p_value)) %>% 
  dplyr::slice(1:5) %>% 
  ungroup() %>% 
  mutate(term_name = tolower(term_name)) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  mutate(term_name = fct_reorder2(term_name, cluster,negative_log10_of_adjusted_p_value)) 

p = ggplot(pathways, aes(x = fct_reorder2(term_name, cluster,negative_log10_of_adjusted_p_value), 
                         y = negative_log10_of_adjusted_p_value, fill = source))
p + geom_bar(stat= "identity") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#31A9B8", "#CF3721", "#258039", "#F5BE41")) +
  xlab("")+
  ylab(expression(paste(-log[10],"(padj)"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        axis.text=element_text(size=8), legend.position = "right", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(col = vars(cluster), scales = "free", space = "free")

