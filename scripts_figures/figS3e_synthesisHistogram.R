rm(list = ls())

library(tidyverse)
library(RColorBrewer)

# upload synthesis rates (crop values above 1.5 for plotting)
# get quartiles
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  mutate(quartile = factor(ntile(cropped, 4)))

# plot histogram
p = ggplot(synthesisRates, aes(x = cropped))
p + geom_histogram(bins = 25,color = "#CF3721", fill = "#CF3721", alpha = 0.9) +
  theme_classic() +
  xlab("synthesis rate [mRNA/cell/min]") +
  ylab("n") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
