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

noTBP = lm(data = plotResTime, formula = synthesis_perCell_perMin ~TFIIA + TFIIB + TFIIF + TFIIE)
summary(noTBP)
noTBPmodel = data.frame(factorName = names(noTBP$coefficients[-1]), p = summary(noTBP)$coefficients[,4][-1], 
                   beta = noTBP$coefficients[-1], type = rep("noTBP", length(noTBP$coefficients[-1]))) 
model = rbind(model, noTBPmodel)

addNaTBP= data.frame(factorName = "TBP", p = NA, 
                     beta = NA, type = "noTBP")
model = rbind(model, addNaTBP)

plotModel = model %>% 
  mutate(negLog10p = -log10(p))

p = ggplot(plotModel, aes(x = factorName, y = beta, fill = type))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("#CF3721", "#31A9B8", "#F5BE41")) +
  xlab("") +
  geom_hline(yintercept = 0)
ggsave("figures/panels/random/lm_beta.pdf", width = 10, height = 6, units = "cm")

p = ggplot(plotModel %>% mutate(negLog10p = ifelse(negLog10p > 15, 15, negLog10p)), aes(x = factorName, y = negLog10p, fill = type))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("#CF3721", "#31A9B8", "#F5BE41")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  xlab("")
ggsave("figures/panels/random/lm_p.pdf", width = 10, height = 6, units = "cm")

write.csv(plotModel, "figures/panels/random/lm.csv", quote = F, row.names = F)

#####################################################################
##### remove the <1 values
#####################################################################
rm(list = ls())

set.seed(42)

# upload synthesis rates for coloring
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  dplyr::filter(resTime != "<1") %>% 
  mutate(resTimesNum = as.numeric(resTime)) %>% 
  dplyr::select(-c(peakName, resTime)) %>% 
  pivot_wider(names_from = "factorName", values_from = "resTimesNum") %>% 
  dplyr::select(gene, TFIIB, TFIIE, TFIIF, TFIIA, TBP) %>% 
  left_join(synthesisRates) 


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

noTBP = lm(data = plotResTime, formula = synthesis_perCell_perMin ~TFIIA + TFIIB + TFIIF + TFIIE)
summary(noTBP)
noTBPmodel = data.frame(factorName = names(noTBP$coefficients[-1]), p = summary(noTBP)$coefficients[,4][-1], 
                        beta = noTBP$coefficients[-1], type = rep("noTBP", length(noTBP$coefficients[-1]))) 
model = rbind(model, noTBPmodel)

addNaTBP= data.frame(factorName = "TBP", p = NA, 
                     beta = NA, type = "noTBP")
model = rbind(model, addNaTBP)

plotModel = model %>% 
  mutate(negLog10p = -log10(p))

p = ggplot(plotModel %>% mutate(beta = ifelse(beta < -0.1, -0.1, beta)), aes(x = factorName, y = beta, fill = type))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("#CF3721", "#31A9B8", "#F5BE41")) +
  xlab("") +
  geom_hline(yintercept = 0)
ggsave("figures/panels/random/lm_beta_noFastSites.pdf", width = 10, height = 6, units = "cm")

p = ggplot(plotModel %>% mutate(negLog10p = ifelse(negLog10p > 15, 15, negLog10p)), aes(x = factorName, y = negLog10p, fill = type))
p + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("#CF3721", "#31A9B8", "#F5BE41")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  xlab("")
ggsave("figures/panels/random/lm_p_noFastSites.pdf", width = 10, height = 6, units = "cm")

write.csv(plotModel, "figures/panels/random/lm_noFastSites.csv", quote = F, row.names = F)
