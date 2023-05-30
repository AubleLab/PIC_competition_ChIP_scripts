rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# upload reasidence times of Rap1
fkh1_rap1 = read.csv("data/Fkh1_Rap1_resTimes.csv") %>% 
  dplyr::rename(Fkh1Rap1 = factorName, resFR = resTime)

# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-c(peakName)) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  mutate(resTimesNumCropped = ifelse(resTimesNum > 20, 20, resTimesNum)) %>% 
  dplyr::select(gene, factorName, resTimesNum, resTimesNumCropped) %>% 
  inner_join(fkh1_rap1) %>% 
  mutate(resFR = ifelse(resFR > 20 & Fkh1Rap1 == "Rap1", 20, resFR))
plotResTime$factorName = factor(plotResTime$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                   "TFIIF", "TFIIE"))

# get correlation coefficients for Rap1 vs factor
corelations = plotResTime %>% 
  group_by(factorName, Fkh1Rap1) %>% 
  summarise(cor = cor(resTimesNum, resFR), 
            rho = cor(resTimesNum, resFR, method = "spearman")) %>% 
  mutate(myLabel = paste0("r ~ ", round(cor, digits = 2)))

# add color for each factor
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(plotResTime, aes(x = resTimesNumCropped, y = resFR, color = factorName))
p + geom_point() +
  theme_bw() +
  #geom_smooth(method = "lm") +
  facet_grid(Fkh1Rap1~factorName, scales = "free_y") +
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
  ylab("Fkh1 / Rap1 residence time [min]") +
  xlab("residence time [min]") +
  geom_text(data = corelations, aes(label = myLabel), y = 20, x = 10, color = "black", size = 3) +
  geom_vline(xintercept = 1, color = "grey40", linetype = "dashed")
ggsave("figures/panels/random/resTime_vs_Rap1Fkh1resTime.pdf", width = 17, height = 9, units = "cm")


# compare Ra1 residence times for long- vs -shortlived TFIIF
longVsShort = plotResTime %>% 
  dplyr::filter(factorName == "TFIIF") %>% 
  mutate(TFIIF_class = ifelse(resTimesNum > 5, "long-lived", "short-lived"))

longVsShort$TFIIF_class = factor(longVsShort$TFIIF_class, levels = c("short-lived", "long-lived"))

p = ggplot(longVsShort, aes(x = TFIIF_class, y = resFR, fill = TFIIF_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.5), 
               aes(color = TFIIF_class), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA, color = "black", alpha = 0.5) +
  theme_classic() +
  xlab("TFIIF class") +
  ylab("Fkh1/Rap1 residence time [min]") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  scale_fill_manual(values= c("#E29930", "#217CA3")) +
  scale_color_manual(values= c("#E29930", "#217CA3")) + 
  stat_compare_means(method = "t.test", label = "p.signif", label.y = 20) +
  facet_wrap(~Fkh1Rap1, scales = "free_y")
ggsave("figures/panels/random/Rap1_FKh1_vs_TFIIFlong_short.pdf", width = 5, height = 6, units = "cm")






