rm(list = ls())

synthesis = read.csv("data/synthesisRatesEstimates_20_vs_60_min.csv") %>% 
  dplyr::filter(!is.na(synthesis60) & !is.na(synthesis20))

r = cor(synthesis$synthesis20, synthesis$synthesis60)
plotSynthesis = synthesis %>% 
  dplyr::select(synthesis20, synthesis60) %>% 
  mutate(cropped20 = ifelse(synthesis20 > 1.5, 1.5, synthesis20)) %>% 
  mutate(cropped60 = ifelse(synthesis60 > 1.5, 1.5, synthesis60))

p = ggplot(plotSynthesis, aes(x = cropped20, y = cropped60))
p + geom_point(alpha = 0.5) +
  coord_fixed() + 
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  xlab(paste("20 min in galactose","synthesis rate [mRNA/cell/min]", sep = "\n"))+
  ylab(paste("60 min in galactose","synthesis rate [mRNA/cell/min]", sep = "\n")) +
  ggtitle(paste0("r = ", round(r, digits = 2)))
ggsave("figures/panels/figS3/synthesis_60_vs20_min.pdf", width = 6, height = 6, units = "cm")
