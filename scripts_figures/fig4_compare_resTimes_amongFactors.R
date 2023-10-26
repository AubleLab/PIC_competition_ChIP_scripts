rm(list = ls())

library(tidyverse)
library(GGally)
library(RColorBrewer)
library(viridis)

set.seed(42)

# upload synthesis rates for coloring
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

combinations = resTimes %>% 
  dplyr::select(gene, factorName) %>% 
  mutate(factorName_x = factor(factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                      "TFIIF", "TFIIE"))) %>% 
  mutate(factorName_y = factor(factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                      "TFIIF", "TFIIE"))) %>% 
  tidyr::expand(gene, factorName_x, factorName_y)

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  mutate(resTimesNumCropped = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  dplyr::select(-c(peakName, resTime, randNum, resTimesNum))

plotData = combinations %>% 
  left_join(plotResTime, by = c("gene" = "gene", "factorName_x" = "factorName")) %>% 
  dplyr::rename(x = resTimesNumCropped) %>% 
  left_join(plotResTime, by = c("gene" = "gene", "factorName_y" = "factorName")) %>% 
  dplyr::rename(y = resTimesNumCropped) %>% 
  left_join(synthesisRates)

plotData$factorName_x = factor(plotData$factorName_x, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                 "TFIIF", "TFIIE"))
plotData$factorName_y = factor(plotData$factorName_y, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                 "TFIIF", "TFIIE"))

p = ggplot(plotData, aes(x = x, y = y))
p + geom_point(alpha = 0.2, color = "grey40") +
  theme_bw() +
  facet_grid(factorName_y ~ factorName_x) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("residence time [min]") +
  ylab("residence time [min]") + coord_fixed()


  
p = ggplot(plotData, aes(x = x, y = y, z = cropped))
p + geom_point(alpha = 0.2, color = "grey70") +
  stat_summary_2d(data = plotData %>% dplyr::filter(!is.na(cropped)), 
                  bins = 20, color = "transparent") +
  theme_bw() +
  facet_grid(factorName_y ~ factorName_x) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  xlab("residence time [min]") +
  ylab("residence time [min]") +
  scale_fill_viridis_c(option = "B", direction = -1)+
  coord_fixed()




