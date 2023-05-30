rm(list = ls())

library(tidyverse)

# load imputed values
imp = as.matrix(read.csv("data/imputed_missing_values/imputed_residence_times.csv") %>% 
                  column_to_rownames(var = "X"))

# load synthesis rates and get quartiles 
TATA = read.csv("data/TATA_information.csv")


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

# a) imputed observations
output = genePCA(inputData = imp, scaled = TRUE)
plotData = output$PCA_out %>% 
  left_join(TATA) 
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
# ggsave("figures/panels/fig4/PCA_genes.pdf", width = 7, height = 7, units = "cm")

p = ggplot(plotData %>% filter(!is.na(TATA_class)), aes(x = PC1, y = PC2, color = TATA_class))
p + geom_point(data = plotData %>% filter(!is.na(TATA_class) & TATA_class == "TATA-less"),alpha = 0.5) +
  geom_point(data = plotData %>% filter(!is.na(TATA_class) & TATA_class == "TATA-containing"),alpha = 0.5) +
  #stat_summary_2d(bins = 30, color = "transparent") +
  theme_bw() +
  xlab(paste('PC1 (', variance[1],'%)', sep = ""))+
  ylab(paste('PC2 (', variance[2],'%)', sep = "")) +
  coord_fixed() +
  scale_color_manual(values= c("#F5BE41", "#31A9B8", "grey70")) +
  ylim(-3, 10)+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=9),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top")
ggsave("figures/panels/TATA_PCA_genes.pdf", width = 6.5, height = 6.5, units = "cm")
