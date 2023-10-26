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
plotData = output$PCA_out %>% 
  left_join(synthesisRates)
variance = output$variance

# p = ggplot(plotData, aes(x = PC1, y = PC2, color = cropped))
# p + geom_point(alpha = 0.1, color = "grey70", size = 0.5) +
#   geom_point(data = plotData %>% dplyr::filter(!is.na(synthesis_perCell_perMin)), 
#              alpha = 0.7, size =0.5) +
#   theme_bw() +
#   xlab(paste('PC1 (', variance[1],'%)', sep = ""))+
#   ylab(paste('PC2 (', variance[2],'%)', sep = "")) +
#   coord_fixed() +
#   scale_color_viridis_c(option = "B", direction = -1) +
#   ylim(-3, 10) +
#   theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent"),
#         text = element_text(size=10),
#         strip.background=element_rect(colour="black",
#                                       fill="white"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top") +
#   guides(color = guide_colorbar(title = paste("synthesis rate", "[mRNA/cell/min]", sep = "\n")))

p = ggplot(plotData, aes(x = PC1, y = PC2, z = cropped))
p + geom_point(alpha = 0.2, color = "grey70") +
  stat_summary_2d(data = plotData %>% dplyr::filter(!is.na(synthesis_perCell_perMin)), 
                  bins = 40, color = "transparent") +
  theme_bw() +
  xlab(paste('PC1 (', variance[1],'%)', sep = ""))+
  ylab(paste('PC2 (', variance[2],'%)', sep = "")) +
  coord_fixed() +
  scale_fill_viridis_c(option = "B", direction = -1)+
  #ylim(-3, 10)+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=9),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top") +
  guides(fill= guide_colorbar(title = paste("mean synthesis rate", "[mRNA/cell/min]", sep = "\n")))
