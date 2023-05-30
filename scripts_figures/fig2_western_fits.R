rm(list = ls())

library(tidyverse)

# upload western blots
western = read.csv("data/westerns.csv")
westernFits = read.csv("data/westerns_fitParam.csv")

western$factorName = factor(western$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))
westernLabels = westernFits %>% 
  mutate(myLabel = paste0(factorName,", n=", round(n_western, digits = 2),  " ,t[1/2]ind=", round(Km_western, digits = 2)))
rownames(westernLabels) = westernLabels$factorName

myLabels =westernLabels[c("TBP", "TFIIA", "TFIIB", "TFIIF", "TFIIE"), "myLabel"]
Wmax = 1
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(western, aes(x = time, y = ratio, color = factorName))
p + geom_smooth(aes(group = factorName),
                method="nls",
                method.args = list(formula = y ~ Wmax*x^n / (Jm^n + x^n), 
                                   start=list(Jm=40, n=3)),
                se=FALSE)  + 
  theme_classic()+ 
  ylab("HA/Myc") +
  scale_color_manual(values = colorPalette, labels = c(bquote("TBP, n=4.95, "~t[1/2]*"ind=33.26"),
                                                       bquote("TFIIA, n=3.85, "~t[1/2]*"ind=57.96"),
                                                       bquote("TFIIB, n=4.28, "~t[1/2]*"ind=42.61"),
                                                       bquote("TFIIF, n=6.58, "~t[1/2]*"ind=42.25"),
                                                       bquote("TFIIE, n=2.94, "~t[1/2]*"ind=40.34")))+
  guides(color = guide_legend(title = "  "))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=10), 
        legend.background = element_rect(fill = "transparent")) +
  geom_hline(yintercept = 0.5) +
  xlab("time [min]")
ggsave("figures/panels/fig2/inductionCurves.pdf", width = 10, height = 4.5, units = "cm")
