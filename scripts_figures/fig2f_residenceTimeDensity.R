rm(list = ls())

library(tidyverse)
set.seed(42)

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(factorName, resTimesNum) %>% 
  mutate(resTimesNum = ifelse(resTimesNum > 20, 20, resTimesNum))

plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                           "TFIIF", "TFIIE"))
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(plotResTime, aes(resTimesNum, color = factorName, fill = factorName))
p + geom_density(alpha = 0.8, lwd = 1) +
  theme_classic() +
  scale_color_manual(values = colorPalette) +
  scale_fill_manual(values = colorPalette) +
  facet_wrap(~factorName, nrow = 1, scales = "free_y") + 
  theme(legend.position = "none") +
  xlab("residence time [min]")+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=10)) +
  geom_vline(xintercept = 1, color = "grey40", linetype = "dashed")
