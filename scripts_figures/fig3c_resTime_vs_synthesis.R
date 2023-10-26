rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# load synthesis rates and get quartiles 
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(quartile = factor(ntile(synthesis_perCell_perMin, 4)))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  inner_join(synthesisRates)
  
plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

p = ggplot(plotResTime, aes(x = factorName, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  ylim(0,20) 

##### ------- add significance testing ------


my_comparisons <- list(c("1", "2"), c("2", "3"), c("3", "4"), 
                       c("1", "3"), c("2", "4"), c("1", "4") )
# 1) TFIIE
TFIIE = plotResTime %>% 
  dplyr::filter(factorName == "TFIIE")
p = ggplot(TFIIE, aes(x = quartile, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  #ylim(0,20)  + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25) 

# 2) TBP
TBP = plotResTime %>% 
  dplyr::filter(factorName == "TBP")
p = ggplot(TBP, aes(x = quartile, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  #ylim(0,20)  + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25) 


# 3) TFIIA
TFIIA = plotResTime %>% 
  dplyr::filter(factorName == "TFIIA")
p = ggplot(TFIIA, aes(x = quartile, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  #ylim(0,20)  + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25)


# 4) TFIIB
TFIIB = plotResTime %>% 
  dplyr::filter(factorName == "TFIIB")
p = ggplot(TFIIB, aes(x = quartile, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  #ylim(0,20)  + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25)


# 5) TFIIF
TFIIF = plotResTime %>% 
  dplyr::filter(factorName == "TFIIF")
p = ggplot(TFIIF, aes(x = quartile, y = resTimesNum, fill = quartile))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = quartile), alpha = 0.5) + 
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
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  #ylim(0,20)  + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25) 

