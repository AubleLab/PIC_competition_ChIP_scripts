rm(list = ls())

library(tidyverse)

set.seed(42)

# load synthesis rates and get quartiles 
synthesisRates = read.csv("data/synthesisRates.csv") 

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  inner_join(synthesisRates) %>% 
  mutate(TR_efficiency = synthesis_perCell_perMin * resTimesNum)

plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(plotResTime, aes(x = factorName, y = log2(TR_efficiency), fill = factorName))
p + geom_violin() +
  theme_classic() +
  xlab(" ") +
  ylab((expression(paste(log[2],"(mRNA per binding event)")))) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #facet_wrap(~factorName, scales = "free", nrow = 1) +
  scale_color_manual(values = colorPalette) +
  scale_fill_manual(values = colorPalette) +
  stat_summary(fun=median, geom="point", shape=16, size=1)

mean_and_median = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  #mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  inner_join(synthesisRates) %>% 
  mutate(TR_efficiency = synthesis_perCell_perMin * resTimesNum) %>% 
  group_by(factorName) %>% 
  summarise(mean = mean(TR_efficiency), med = median(TR_efficiency))

