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
  mutate(resTimesNumCropped = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  inner_join(synthesisRates) %>% 
  mutate(croppedRate = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))
plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

correlations = plotResTime %>% 
  group_by(factorName) %>% 
  summarise(correlation = cor(synthesis_perCell_perMin, resTimesNum)) %>% 
  mutate(myLabel = paste0("r ~ ", round(correlation, digits = 2)))
correlations$factorName = factor(correlations$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(plotResTime, aes(x = resTimesNumCropped, y = croppedRate, color = factorName))
p + geom_point(alpha = 0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + 
  theme_classic() +
  ylab(paste("synthesis rate", "[mRNA/cell/min]", sep = "\n")) +
  xlab("residence time [min]") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(data = correlations, aes(label = myLabel), y = 1.7, x = 10, color = "black", size = 3) +
  facet_wrap(~factorName, nrow = 1) +
  scale_color_manual(values = colorPalette) +
  #geom_smooth(method = "lm", color = "black") +
  ylim(0, 1.7)
ggsave("figures/panels/figSc/resTime_vs_synthesisRate.pdf", width = 17, height = 4.5, units = "cm")

