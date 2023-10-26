rm(list = ls())

library(tidyverse)
library(tidytext)
library(svglite)

set.seed(42)

# load imputed values
imp = as.matrix(read.csv("data/imputed_missing_values/imputed_residence_times_fastSitesExcluded.csv") %>% 
                  column_to_rownames(var = "X"))

# upload synthesis rates
synthesis = read.csv("data/synthesisRates.csv")



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
  inner_join(synthesis) %>% 
  pivot_longer(cols = PC1:PC5, names_to = "PC", values_to = "PC_value") %>% 
  group_by(PC) %>% 
  summarise(corVal = cor(PC_value, synthesis_perCell_perMin)) %>% 
  drop_na()

p = ggplot(corData %>% dplyr::filter(PC %in% c("PC1", "PC2")), aes(x = PC, y = corVal, fill = corVal))
p + geom_bar(stat = "identity", color = "grey70", lwd = 0.3) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=10), legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #guides(fill="none", color = "none") +
  ylab("r") +
  xlab("")  +
  scale_fill_gradientn(colors = c("#32384D", "white","#E29930"), 
                       limits= c(-1,1)) +
  ylim(-1, 1)
