rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ggrepel)


# upload table with residence times and replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(c(gene, TFIIF)) %>% 
  drop_na() %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(TFIIF = ifelse(TFIIF == "<1", randNum, as.numeric(TFIIF))) %>% 
  dplyr::select(gene, TFIIF)


# upload residence time table - select TFIIF 
# and separate those into long-lived (5 min thresh.) or not
# extract the short-lived TFIIF genes
shortTFIIF = plotResTime %>% 
  mutate(shortTFIIF = ifelse(TFIIF <5, "yes", "no")) %>% 
  mutate(shortTFIIF = factor(shortTFIIF, levels = c("yes", "no")))
#write.table(shortTFIIF$gene,"data/analysis/short_TFIIF/genes.txt", quote = F,row.names = F, col.names = F)

# paste the lists of long lived and short lived genes to g:Profiler
# upload results from g:Profiler: TF targets (for duplicated results from g:Profiler 
# pick the one with the most significant enrichment)
gprofiler = read.csv("data/analysis/short_TFIIF/gProfiler.csv") %>% 
  dplyr::filter(source == "TF") %>% 
  separate(term_name, into = c("toss", "TF"), sep = ": ", extra = "drop") %>% 
  dplyr::select(-toss) %>% 
  separate(TF, into = c("TF", "toss"), sep = ";", extra = "drop") %>% 
  dplyr::select(-toss) %>% 
  mutate(TF = ifelse(endsWith(TF, "p"), str_sub(TF, end = -2), TF)) %>% 
  group_by(TF) %>% 
  dplyr::slice(which.max(negative_log10_of_adjusted_p_value)) %>% 
  dplyr::select(TF, negative_log10_of_adjusted_p_value)%>% 
  mutate(TF = str_to_title(TF))

# compare results from  significant TF target enrichment (Yeats Epigeneome database gene targets): 
# these were generated analogously with "figS7f_TFIIF_shortLived_YeastEpi_enrichment.R") 
# and compare to results from g:Profiler
sigRes = read.csv("data/analysis/short_TFIIF/sigTFenrichment.csv")
plotRes = sigRes %>% 
  dplyr::select(TF, logFDR) %>% 
  mutate(TF = str_to_title(TF)) %>% 
  full_join(gprofiler) %>% 
  replace_na(list(negative_log10_of_adjusted_p_value = 0, logFDR = 0))

# plot results from g:Profiler "validated" with our TF enrichment results
p = ggplot(plotRes, aes(x = logFDR, y = negative_log10_of_adjusted_p_value))
p +geom_point() +
  theme_classic() +
  geom_label_repel(aes(label = TF),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50', size = 2.5) +
  ylab(expression(paste("TRANSFAC: ",-log[10],"(padj)")))+
  xlab(expression(paste("Yeast Epi.: ",-log[10],"(FDR)"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9))


# upload pathway results from g:Profiler (get maximum top 5 from each database)
pathways = read.csv("data/analysis/short_TFIIF/gProfiler.csv") %>% 
  dplyr::filter(source %in% c("GO:BP", "GO:MF", "KEGG", "WP")) %>% 
  group_by(source) %>% 
  arrange(desc(negative_log10_of_adjusted_p_value)) %>% 
  dplyr::slice(1:5) %>% 
  ungroup() %>% 
  mutate(term_name = tolower(term_name)) %>% 
  mutate(term_name = fct_reorder(term_name,negative_log10_of_adjusted_p_value)) 

#plot
p = ggplot(pathways, aes(y = fct_reorder(term_name,negative_log10_of_adjusted_p_value), 
                         x = negative_log10_of_adjusted_p_value, fill = source))
p + geom_bar(stat= "identity") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#31A9B8", "#CF3721", "#258039", "#F5BE41")) +
  ylab("")+
  xlab(expression(paste(-log[10],"(padj)"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        axis.text=element_text(size=8), legend.position = "top") 
