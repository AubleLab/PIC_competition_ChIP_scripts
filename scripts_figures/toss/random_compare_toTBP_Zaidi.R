rm(list = ls())

library(tidyverse)
library(viridis)

set.seed(42)
# upload the published data
Zaidi = read.csv("/Users/lilcrusher/competitionChIP/TBP_Zaidi/TBPresTime_Zaidi_etal_SciRep_2017.csv")

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(peakName:TBP) %>% 
  dplyr::filter(!is.na(TBP))

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(TBPNum = ifelse(TBP == "<1", randNum, as.numeric(TBP))) %>% 
  select(-randNum) %>% 
  inner_join(Zaidi, b = c("gene" = "Locus")) %>% 
  select(-X) %>% 
  dplyr::rename(Zaidi = resTime)


r = cor(plotResTime$TBPNum, plotResTime$Zaidi)
# plot
p = ggplot(plotResTime, aes(y = TBPNum, x = Zaidi))
p + geom_bin2d(binwidth = c(.3, .3)) +
  theme_bw() +
  coord_fixed() +
  ylab(paste("TBP res. time [min]","this study", sep = "\n")) +
  xlab(paste("TBP res. time [min]", "Zaidi et al, 2017", sep = "\n")) +
  scale_fill_viridis() +
  ylim(0,5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))
ggsave("figures/panels/random/TBP_vs_Zaidi_density.pdf", width = 5.8, height = 2.5)

p + geom_point(alpha = 0.4, color = "grey60") +
  theme_bw() +
  coord_fixed() +
  ylab(paste("TBP res. time [min]","this study", sep = "\n")) +
  xlab(paste("TBP res. time [min]", "Zaidi et al, 2017", sep = "\n")) +
  ylim(0,5) +
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))
ggsave("figures/panels/random/TBP_vs_Zaidi.pdf", width = 5, height = 2)
