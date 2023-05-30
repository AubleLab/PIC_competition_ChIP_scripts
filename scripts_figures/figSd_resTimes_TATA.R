rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# load synthesis rates and get quartiles 
TATA = read.csv("data/TATA_information.csv") 
  

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  inner_join(TATA)

plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))
plotResTime$TATA_class = factor(plotResTime$TATA_class, levels = c("TATA-less", "TATA-containing"))

p = ggplot(plotResTime, aes(x = factorName, y = resTimesNum, fill = TATA_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = TATA_class), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA, color = "black") +
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
  scale_fill_manual(values= c("#F5BE41", "#31A9B8")) +
  scale_color_manual(values= c("#F5BE41", "#31A9B8")) +
  ylim(0,20) + stat_compare_means(method = "t.test", label = "p.signif", label.y = 20)
ggsave("figures/panels/figSd/resTime_vs_TATA.pdf", width = 6, height = 6, units = "cm")


