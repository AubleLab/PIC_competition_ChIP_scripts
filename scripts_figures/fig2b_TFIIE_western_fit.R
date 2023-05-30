rm(list = ls())

library(tidyverse)

# upload western blots
western = read.csv("data/westerns.csv") %>% 
  dplyr::filter(factorName == "TFIIE")
westernFits = read.csv("data/westerns_fitParam.csv") %>% 
  dplyr::filter(factorName == "TFIIE")

westernLabels = westernFits %>% 
  mutate(myLabel = paste0(factorName,", n=", round(n_western, digits = 2),  " ,t[1/2]ind=", 
                          round(Km_western, digits = 2)))
rownames(westernLabels) = westernLabels$factorName

Wmax = 1

p = ggplot(western, aes(x = time, y = ratio))
p + geom_smooth(method="nls",
                method.args = list(formula = y ~ Wmax*x^n / (Jm^n + x^n), 
                                   start=list(Jm=40, n=3)),
                se=FALSE, lwd = 0.7, color = "#E64A45")  + 
  theme_classic()+ 
  ylab("HA/Myc") +
  xlab("time [min]") +
  annotate("text", y= 0.05, x = 75, label = deparse(bquote("n=2.94, "~t[1/2]*"ind=40.34")), 
           parse = TRUE, color = "black", size = 3) +
  guides(color = guide_legend(title = "  "))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=10), 
        legend.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size = 9)) +
  geom_point(size = 0.8, color = "#E64A45") +
  geom_errorbar(aes(ymin=ratio-ratio_sd, ymax=ratio+ratio_sd), 
                position=position_dodge(0.05), color = "#E64A45") +
  ggtitle("TFIIE")
#ggsave("figures/panels/fig2/TFIIE_inductionCurve.pdf", width = 5.2, height = 4.7, units = "cm")
