rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# get list of TATA promoters
TATA = read.csv("data/TATA_information.csv") 


# load synthesis rates and attach info about TATA promoters
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  inner_join(TATA) %>% 
  mutate(x= "X")


synthesisRates$TATA_class = factor(synthesisRates$TATA_class, levels = c("TATA-less", "TATA-containing"))

# plot
p = ggplot(synthesisRates, aes(x = x , y = cropped, fill = TATA_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = TATA_class), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA, color = "black") +
  theme_classic() +
  xlab(" ") +
  ylab(paste("synthesis rate", "[mRNA/cell/min]", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #facet_wrap(~factorName, scales = "free", nrow = 1) +
  scale_fill_manual(values= c("#F5BE41", "#31A9B8")) +
  scale_color_manual(values= c("#F5BE41", "#31A9B8")) +
  ylim(0,1.6) + 
  stat_compare_means(label.y = 1.6)


