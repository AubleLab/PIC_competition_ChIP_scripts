rm(list = ls())


library(tidyverse)
library(RColorBrewer)

# upload published synthesis rates
TR = read.delim("data/previously_published/Garcia_transcription_rates_galactose.txt") %>% 
  dplyr::select(-gene_name) %>% 
  dplyr::rename(gene = ORF)%>% 
  mutate(cropped_TR = ifelse(TR_galactose > .5, .5, TR_galactose))

# upload synthesis rates (crop values above 1.5 for plotting)
# get quartiles
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin)) %>% 
  inner_join(TR)

rho = cor(synthesisRates$synthesis_perCell_perMin, 
        synthesisRates$TR_galactose, use = "complete.obs", 
        method = "spearman")

r = cor(synthesisRates$synthesis_perCell_perMin, 
        synthesisRates$TR_galactose, use = "complete.obs", 
        method = "pearson")

p = ggplot(synthesisRates, aes(y = cropped_TR, x = cropped))
p + geom_bin2d(bins = 50) +
  theme_classic() +
  #coord_fixed() +
  xlab(paste("synthesis rate [mRNA/cell/min]","this study", sep = "\n")) +
  ylab(paste("TR [mol/min] exp. YPGAL", "García-Martínez et al, 2004", sep = "\n")) +
  scale_fill_viridis() +
  theme(panel.background = element_rect(fill = "transparent", color = "black"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), legend.position = "top") 



