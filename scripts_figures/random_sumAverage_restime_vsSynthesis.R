rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)


set.seed(42)


# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-c(peakName, TBP)) %>% 
  pivot_longer(cols = TFIIA:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum) %>% 
  drop_na() %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  drop_na() %>% 
  pivot_longer(cols = TFIIA:TFIIF, names_to = "factorName", values_to = "resTime") %>% 
  group_by(gene) %>% 
  summarise(meanRes = mean(resTime), sumRes = sum(resTime))

# upload synthesis rates
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin)) %>% 
  inner_join(plotResTime) %>% 
  distinct() %>% 
  pivot_longer(cols = meanRes:sumRes, names_to = "statName", values_to = "value")

r = synthesisRates %>% 
  group_by(statName) %>% 
  summarise(r = cor(value, synthesis_perCell_perMin), rho = cor(value, synthesis_perCell_perMin, method = "spearman"))


p = ggplot(synthesisRates, aes(x = value, y = cropped))
p + geom_point(alpha = 0.3) + 
  theme_classic() +
  geom_text(data = r, aes(label = paste0("r = ",round(r, digits = 2))), x = 10, y = 1.7) +
  geom_text(data = r, aes(label = paste0("rho = ",round(rho, digits = 2))), x = 10, y = 1.6) +
  facet_wrap(~statName, scales = "free") +
  ylim(0, 1.8) +
  xlab("mean or sum of residence time")
ggsave("figures/panels/random/mean_and_sum_resTime_vs_synthesis.pdf", width = 7, height = 4)



