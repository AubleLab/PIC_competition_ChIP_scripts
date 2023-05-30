rm(list = ls())

library(tidyverse)
library(RColorBrewer)

# upload synthesis rates (crop values above 1.5 for plotting)
# get quartiles
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  mutate(quartile = factor(ntile(cropped, 4)))

# plot density plot of quartiles
p = ggplot(synthesisRates, aes(x = cropped, fill = quartile))
p + geom_density(color = "black", alpha = 0.9) +
  theme_classic() +
  xlab("synthesis rate [mRNA/cell/min]") +
  ylab("density") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Reds")
#ggsave("figures/panels/fig3/synthesisRates_quartiles.pdf", width = 5, height = 4, units = "cm")

