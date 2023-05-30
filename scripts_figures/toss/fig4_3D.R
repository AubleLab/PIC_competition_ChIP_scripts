rm(list = ls())

library(tidyverse)
library(GGally)
library(RColorBrewer)
library(viridis)
library(gg3D)
library(cowplot)

set.seed(42)

# upload synthesis rates for coloring
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  mutate(resTimesNumCropped = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  dplyr::select(-c(peakName, resTime, randNum, resTimesNum)) %>% 
  pivot_wider(names_from = "factorName", values_from = "resTimesNumCropped") %>% 
  dplyr::select(gene, TFIIB, TFIIE, TFIIF) %>% 
  left_join(synthesisRates) %>% 
  drop_na()

p = ggplot(plotResTime, aes(y=TFIIB, x=TFIIF, z=TFIIE, color=cropped))
p + axes_3D() +
  stat_3D()+ 
  scale_color_viridis_c(option = "B", direction = -1) +
  labs_3D(
    labs=c("TFIIF", "TFIIB", "TFIIE"),
    hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90)) +
  theme_void() +
  coord_fixed()
ggsave("figures/panels/fig4/3D_1.pdf", width = 10, height = 8, units = "cm")

p = ggplot(plotResTime, aes(z=TFIIB, x=TFIIF, y=TFIIE, color=cropped))
p + axes_3D() +
  stat_3D()+ 
  scale_color_viridis_c(option = "B", direction = -1) +
  labs_3D(
    labs=c("TFIIF", "TFIIE", "TFIIB"),
    hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90)) +
  theme_void() +
  coord_fixed()
ggsave("figures/panels/fig4/3D_2.pdf", width = 10, height = 8, units = "cm")


p = ggplot(plotResTime, aes(x=TFIIB, z=TFIIF, y=TFIIE, color=cropped))
p + axes_3D() +
  stat_3D()+ 
  scale_color_viridis_c(option = "B", direction = -1) +
  labs_3D(
    labs=c("TFIIB", "TFIIE", "TFIIF"),
    hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90)) +
  theme_void() +
  coord_fixed()
ggsave("figures/panels/fig4/3D_3.pdf", width = 10, height = 8, units = "cm")


