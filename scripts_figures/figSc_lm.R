rm(list = ls())

library(tidyverse)
library(GGally)
library(RColorBrewer)
library(viridis)
library(gg3D)
library(cowplot)

set.seed(42)

# upload synthesis rates for coloring
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(-c(peakName, resTime, randNum)) %>% 
  pivot_wider(names_from = "factorName", values_from = "resTimesNum") %>% 
  dplyr::select(gene, TFIIB, TFIIE, TFIIF, TFIIA, TBP) %>% 
  left_join(synthesisRates) 

# get linear models for individual TFs
TFIIE = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TFIIE)
summary(TFIIE)  
model = data.frame(factorName = "TFIIE", p = summary(TFIIE)$coefficients[,4][2], beta = TFIIE$coefficients[2], type = "single")


TFIIF = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TFIIF)
summary(TFIIF) 
TFIIFmodel = data.frame(factorName = "TFIIF", p = summary(TFIIF)$coefficients[,4][2], beta = TFIIF$coefficients[2], type = "single")
model = rbind(model, TFIIFmodel)


TFIIB = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TFIIB)
summary(TFIIB)  
TFIIBmodel = data.frame(factorName = "TFIIB", p = summary(TFIIB)$coefficients[,4][2], beta = TFIIB$coefficients[2], type = "single")
model = rbind(model, TFIIBmodel)


TFIIA = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TFIIA)
summary(TFIIA)  
TFIIAmodel = data.frame(factorName = "TFIIA", p = summary(TFIIA)$coefficients[,4][2], beta = TFIIA$coefficients[2], type = "single")
model = rbind(model, TFIIAmodel)

TBP = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TBP)
summary(TBP)  
TBPmodel = data.frame(factorName = "TBP", p = summary(TBP)$coefficients[,4][2], beta = TBP$coefficients[2], type = "single")
model = rbind(model, TBPmodel)

all = lm(data = plotResTime, formula = synthesis_perCell_perMin ~ TBP + TFIIA + TFIIB + TFIIF + TFIIE)
summary(all)
allModel = data.frame(factorName = names(all$coefficients[-1]), p = summary(all)$coefficients[,4][-1], beta = all$coefficients[-1], 
                      type = rep("combined", length(all$coefficients[-1]))) 
model = rbind(model, allModel)


plotModel = model %>% 
  mutate(negLog10p = -log10(p))

plotModel$factorName = factor(plotModel$factorName, levels = c("TBP", "TFIIA", "TFIIB", "TFIIF", "TFIIE"))
colorPalette = c("#F2C249", "#E6772E", "#3D4C53", "#4DB3B3", "#E64A45")
p = ggplot(plotModel %>% dplyr::filter(type == "single"), aes(x = factorName, y = beta, fill = factorName))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        legend.background = element_rect(fill = "transparent"),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colorPalette) +
  xlab("") +
  geom_hline(yintercept = 0) +
  ylim(-0.1, 0.07)
ggsave("figures/panels/figSc/lm_beta_singleParam.pdf", width = 6, height = 4, units = "cm")


p = ggplot(plotModel %>% dplyr::filter(type == "combined"), aes(x = factorName, y = beta, fill = factorName))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        legend.background = element_rect(fill = "transparent"),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colorPalette) +
  xlab("") +
  geom_hline(yintercept = 0)+
  ylim(-0.1, 0.07)
ggsave("figures/panels/figSc/lm_beta_combinedParam.pdf", width = 6, height = 4, units = "cm")
