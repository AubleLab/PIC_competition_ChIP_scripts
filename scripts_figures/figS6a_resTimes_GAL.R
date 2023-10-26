rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# get GAL genes
GAL = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(locus=V4, gene = V6)%>% 
  dplyr::filter(startsWith(gene, "GAL")) %>% 
  distinct() %>% 
  dplyr::rename(gene=locus, sysName = gene)


# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  left_join(GAL) %>% 
  mutate(GAL_class = ifelse(is.na(sysName), "no", "yes"))

plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))
plotResTime$GAL_class = factor(plotResTime$GAL_class, levels = c("no", "yes"))

p = ggplot(plotResTime, aes(x = factorName, y = resTimesNum, fill = GAL_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = GAL_class)) + 
  geom_boxplot(outlier.shape=NA, color = "black", alpha = 0.5) +
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
  scale_fill_manual(values= c("#E6772E", "#3D4C53")) +
  scale_color_manual(values= c("#E6772E", "#3D4C53")) +
  ylim(0,20) + stat_compare_means(method = "wilcox.test", label = "p.signif", label.y = 20)


