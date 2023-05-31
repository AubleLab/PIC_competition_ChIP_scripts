rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# upload matching systematic names from www.yeastgenome.org:
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)

# upload the residence time table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

# identify ribosomal subunits
plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  left_join(geneConvertTable) %>% 
  mutate(ribo = ifelse(str_detect(gene, "RPL") | str_detect(geneSys, "RPL"), "yes", "no")) %>% 
  replace_na(list(ribo = "no")) %>% 
  distinct()

plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))
plotResTime$ribo = factor(plotResTime$ribo, levels = c("no", "yes"))

p = ggplot(plotResTime, aes(x = factorName, y = resTimesNum, fill = ribo))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = ribo), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA) +
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
  #facet_wrap(~factorName, scales = "free", nrow = 1) +
  scale_fill_manual(values= c("#4DB3B3", "#E64A45")) +
  scale_color_manual(values= c("#4DB3B3", "#E64A45")) +
  ylim(0,20) + stat_compare_means(method = "t.test", label = "p.signif", label.y = 20)
#ggsave("figures/panels/figS6/resTime_vs_RIBO.pdf", width = 6, height = 6, units = "cm")


