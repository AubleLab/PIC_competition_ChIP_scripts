rm(list = ls())

library(tidyverse)

# load imputed values
imp = as.matrix(read.csv("data/imputed_missing_values/imputed_residence_times.csv") %>% 
                  column_to_rownames(var = "X"))

# load synthesis rates and get quartiles 
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  mutate(quartile = factor(ntile(cropped, 4)))


# ======= 1) pca of factors ===========
#2) PCA of genes
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

# a) imputed observations
output = genePCA(inputData = imp, scaled = TRUE)
plotData = output$PCA_out

combinations = plotData %>% 
  dplyr::select(gene:PC5) %>% 
  pivot_longer(cols = PC1:PC5, names_to = "PC", values_to = "value") %>% 
  dplyr::select(gene, PC) %>% 
  mutate(PC_x = factor(PC)) %>% 
  mutate(PC_y = factor(PC)) %>% 
  tidyr::expand(gene, PC_x, PC_y)


plotPC = plotData %>% 
  pivot_longer(cols = PC1:PC5, names_to = "PC", values_to = "value")

plotAll = combinations %>% 
  left_join(plotPC, by = c("gene" = "gene", "PC_x" = "PC")) %>% 
  dplyr::rename(x = value) %>% 
  left_join(plotPC, by = c("gene" = "gene", "PC_y" = "PC")) %>% 
  dplyr::rename(y = value) %>% 
  left_join(synthesisRates)

plotAll$PC_x = factor(plotAll$PC_x, levels = c("PC1", "PC2", "PC3", "PC4", "PC5"))
plotAll$PC_y = factor(plotAll$PC_y, levels = c("PC1", "PC2", "PC3", "PC4", "PC5"))

p = ggplot(plotAll, aes(x = x, y = y))
p + geom_point(alpha = 0.2, color = "grey40") +
  theme_bw() +
  facet_grid(PC_y ~ PC_x) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") + coord_fixed() +
  xlim(-10, 10) +
  ylim(-10, 10)
ggsave("figures/panels/figSd/compare_allPCs.pdf", width = 14, height = 14, units = "cm")



p = ggplot(plotAll, aes(x = x, y = y, z = cropped))
p + geom_point(alpha = 0.2, color = "grey70") +
  stat_summary_2d(data = plotAll %>% dplyr::filter(!is.na(cropped)), 
                  bins = 30, color = "transparent") +
  theme_bw() +
  facet_grid(PC_y ~ PC_x) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  xlab("") +
  ylab("") +
  scale_fill_viridis_c(option = "B", direction = -1)+
  coord_fixed()+
  xlim(-10, 10) +
  ylim(-10, 10)
ggsave("figures/panels/figSd/compare_allPCs_synthesisRate.pdf", width = 14, height = 14, units = "cm")








