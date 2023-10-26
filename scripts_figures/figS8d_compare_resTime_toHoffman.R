rm(list = ls())

library(tidyverse)
library(ggrepel)
set.seed(42)

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum) 

# upload Hoffman data
Hoffman = read.csv("data/previously_published/Hoffman_TBP_estimates.csv") %>% 
  dplyr::select(-X) %>% 
  mutate(Hoffman = as.numeric(resTime)) %>% 
  select(-resTime)


# matching systematic names from www.yeastgenome.org:
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab
geneConvertTable = read.delim("/Users/lilcrusher/annotations/yeast/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)

compareResTime = plotResTime %>% 
  left_join(geneConvertTable) %>% 
  distinct() %>% 
  inner_join(Hoffman)
# 


p = ggplot(compareResTime %>% dplyr::filter(method == "CLKv2"), 
           aes(x = resTimesNum, y = Hoffman, color = factorName, fill = factorName))
p + geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") + 
  geom_point(size = 2) +
  theme_classic() + 
  geom_label_repel(aes(label = geneSys),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   color = "black", 
                   segment.color = 'grey50', size = 2.5) +
  #coord_fixed() +
  xlab(paste("residence time [min]", "this study", sep = "\n")) +
  ylab(paste("residence time [min]", "Zaidi, Hoffman et al, 2017", sep = "\n")) +
  scale_color_manual(values = c("#F2C249", "#E64A45")) +
  scale_fill_manual(values = c("#F2C249", "#E64A45")) +
  xlim(0, 18) +
  ylim(0,6) +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = NA),
        text = element_text(size=9), 
        legend.position = "none")


