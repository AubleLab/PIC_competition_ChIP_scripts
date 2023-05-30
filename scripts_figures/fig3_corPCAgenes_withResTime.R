rm(list = ls())

library(tidyverse)
library(tidytext)
library(svglite)

set.seed(42)

# load imputed values
imp = as.matrix(read.csv("data/imputed_missing_values/imputed_residence_times.csv") %>% 
                  column_to_rownames(var = "X"))

# upload table replace <1 min with a random number between 0-1
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) 

# PCA of genes
genePCA = function(inputData, scaled = TRUE){
  if (scaled == TRUE){
    pca = prcomp(inputData, scale. = TRUE)
  } else {
    pca = prcomp(inputData, scale. = FALSE)
  }
  
  variance = round((pca$sdev^2 / sum(pca$sdev^2))*100, digits = 2)
  PCA_out = as.data.frame(pca$x) %>% 
    rownames_to_column(var = "gene")
  
  output = list(PCA_out = PCA_out, variance = variance)
  
  return(output)
}

# a) PCA ofimputed observations
output = genePCA(inputData = imp, scaled = TRUE)
output$variance


corData = output$PCA_out %>% 
  left_join(plotResTime) %>% 
  dplyr::select(-c(peakName, resTime, randNum)) %>% 
  pivot_longer(cols = PC1:PC5, names_to = "PC", values_to = "PC_value") %>% 
  group_by(PC, factorName) %>% 
  summarise(corVal = cor(PC_value, resTimesNum)) %>% 
  drop_na()

p = ggplot(corData, aes(x = reorder_within(factorName, corVal, PC), y = corVal, fill = corVal))
p + geom_bar(stat = "identity") +
  facet_wrap(~PC, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size = 7)) + 
  #guides(fill="none", color = "none") +
  ylab("r") +
  xlab("")  +
  #scale_fill_gradient2(low = "#32384D", high = "#E29930", mid = "white", midpoint = 0) +
  scale_fill_gradientn(colors = c("#32384D", "white","#E29930"), 
                       limits= c(-1,1)) +
  ylim(-1, 1)
ggsave("figures/panels/fig3/cor_PC_genes_with_resTime.pdf", width = 15, height = 5, units = "cm")
