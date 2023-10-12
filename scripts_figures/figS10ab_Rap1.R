rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# upload residence times of Rap1
rap1 = read.csv("data/previously_published/Rap1_resTime_Lickwar.csv") %>% 
  mutate(gene = ifelse(startsWith(gene,"i"), sub(".", "", gene), gene)) %>% 
  mutate(Rap1 = resTime) %>% 
  mutate(Rap1_z = Z_resTime) %>% 
  dplyr::select(gene, Rap1, Rap1_z)

# upload residence time table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-c(peakName)) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

# attach Rap1 table
plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  mutate(resTimesNumCropped = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  dplyr::select(gene, factorName, resTimesNum, resTimesNumCropped) %>% 
  inner_join(rap1)
plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

# get correlation coefficients for Rap1 vs factor
corelations = plotResTime %>% 
  group_by(factorName) %>% 
  summarise(cor = cor(resTimesNum, Rap1), 
            corZ = cor(resTimesNum, Rap1_z), 
            rho = cor(resTimesNum, Rap1, method = "spearman"), 
            rhoZ = cor(resTimesNum, Rap1_z, method = "spearman")) %>% 
  mutate(myLabel = paste0("r ~ ", round(cor, digits = 2)))

# add color for each factor
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
# plot
p = ggplot(plotResTime, aes(x = resTimesNumCropped, y = Rap1, color = factorName))
p + geom_point() +
  theme_bw() +
  #geom_smooth(method = "lm") +
  facet_wrap(~factorName, nrow = 1) +
  scale_color_manual(values = colorPalette) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        legend.background = element_rect(fill = "transparent"),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab(paste("Rap1 residence time [min]", "Lickwar et al, 2012", sep = "\n")) +
  xlab("residence time [min]") +
  geom_text(data = corelations, aes(label = myLabel), y = 150, x = 10, color = "black", size = 3) +
  geom_vline(xintercept = 1, color = "grey40", linetype = "dashed") +
  ylim(25, 155)
#ggsave("figures/panels/figS8/resTime_vs_Rap1resTime.pdf", width = 17, height = 4.5, units = "cm")


# compare Ra1 residence times for long- vs -shortlived TFIIF
longVsShort = plotResTime %>% 
  dplyr::filter(factorName == "TFIIF") %>% 
  mutate(TFIIF_class = ifelse(resTimesNum > 5, "long-lived", "short-lived"))

longVsShort$TFIIF_class = factor(longVsShort$TFIIF_class, levels = c("short-lived", "long-lived"))
# plot
p = ggplot(longVsShort, aes(x = TFIIF_class, y = Rap1, fill = TFIIF_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.5), 
               aes(color = TFIIF_class), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA, color = "black", alpha = 0.5) +
  theme_classic() +
  xlab("TFIIF class") +
  ylab(paste("Rap1 residence time [min]", "Lickwar et al, 2012", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  scale_fill_manual(values= c("#E29930", "#217CA3")) +
  scale_color_manual(values= c("#E29930", "#217CA3")) + 
  stat_compare_means(label.y = 150)
#ggsave("figures/panels/figS8/Rap1_vs_TFIIFlong_short.pdf", width = 2.8, height = 6, units = "cm")






