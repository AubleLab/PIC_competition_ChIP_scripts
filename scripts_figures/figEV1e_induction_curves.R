rm(list = ls())

library(tidyverse)

# upload 0-1 normalized HA/Myc ratios as measured by western blotting
# along with Hill fit parameters 
western = read.csv("data/westerns.csv")
westernFits = read.csv("data/westerns_fitParam.csv")

# relevel labels for plotting
western$factorName = factor(western$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                           "TFIIF", "TFIIE"))
# make labels for legend
westernLabels = westernFits %>% 
  mutate(myLabel = paste0(factorName,", n=", round(n_western, digits = 2),  " ,t[1/2]ind=", round(Km_western, digits = 2)))
rownames(westernLabels) = westernLabels$factorName

myLabels =westernLabels[c("TBP", "TFIIA", "TFIIB", "TFIIF", "TFIIE"), "myLabel"]
# hard set maximum saturation to 1
Wmax = 1
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
#plot
p = ggplot(western, aes(x = time, y = ratio, color = factorName))
p + geom_smooth(aes(group = factorName),
                method="nls",
                method.args = list(formula = y ~ Wmax*x^n / (Jm^n + x^n), 
                                   start=list(Jm=40, n=3)),
                se=FALSE, lwd = 0.7)  + 
  theme_bw()+ 
  ylab("HA/Myc") +
  xlab("time [min]") +
  scale_color_manual(values = colorPalette, labels = c(bquote("n=4.95, "~t[1/2]*"ind=33.26"),
                                                       bquote("n=3.85, "~t[1/2]*"ind=57.96"),
                                                       bquote("n=4.28, "~t[1/2]*"ind=42.61"),
                                                       bquote("n=6.58, "~t[1/2]*"ind=42.25"),
                                                       bquote("n=2.94, "~t[1/2]*"ind=40.34")))+
  guides(color = guide_legend(title = "  "))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        legend.background = element_rect(fill = "transparent"),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin=ratio-ratio_sd, ymax=ratio+ratio_sd), 
                position=position_dodge(0.05)) +
  facet_wrap(~factorName, nrow = 1)
