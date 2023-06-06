rm(list = ls())

library(tidyverse)
library(viridis)

set.seed(42)
# upload the published data with TBP estimates
Zaidi = read.csv("data/previously_published/TBPresTime_Zaidi_etal_SciRep_2017.csv") %>% 
  dplyr::select(-X)

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_TBPonly_with_tRNA.csv") 

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(TBPNum = ifelse(time == "<1", randNum, as.numeric(time))) %>% 
  select(-randNum) %>% 
  inner_join(Zaidi, b = c("gene" = "Locus")) %>% 
  dplyr::rename(Zaidi = resTime) %>% 
  mutate(geneType = ifelse(startsWith(gene, "t"), "tRNA", 
                           ifelse(startsWith(gene, "s"), "snRNA", "mRNA")))

rFull = cor(plotResTime$TBPNum, plotResTime$Zaidi)
r = plotResTime %>% 
  group_by(geneType) %>% 
  summarise(r = cor(TBPNum, Zaidi))
# plot
colorPalette =  c("#E29930", "#217CA3", "#32384D")

# plot
p = ggplot(plotResTime, aes(y = TBPNum, x = Zaidi, color = geneType))
p + geom_point(alpha = 0.5, size = 1) +
  #  geom_density_2d(color = "grey50", binwidth = 0.015, alpha = 0.7) +
  geom_point(data = plotResTime %>% dplyr::filter(geneType != "mRNA"), size = 1) +
  theme_classic() +
  #coord_fixed() +
  ylab(paste("TBP res. time [min]","this study", sep = "\n")) +
  xlab(paste("TBP res. time [min]", "Zaidi et al, 2017", sep = "\n")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  #facet_wrap(~geneType, nrow = 1)+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = NA),
        text = element_text(size=9), 
        legend.position = "none") +
  scale_color_manual(values = colorPalette)

#ggsave("figures/panels/figS10/TBP_vs_Zaidi_density.pdf", width = 5.5, height = 5.5, units = "cm")

