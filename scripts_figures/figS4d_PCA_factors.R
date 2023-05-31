rm(list = ls())

library(tidyverse)

# load imputed values
imp = as.matrix(read.csv("data/imputed_missing_values/imputed_residence_times.csv") %>% 
                  column_to_rownames(var = "X"))

# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-peakName) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

meanResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(factorName, resTimesNum) %>% 
  group_by(factorName) %>% 
  summarise(meanResTime = mean(resTimesNum))
                
# ======= 1) pca of factors ===========
factorPCA = function(imputedData){
  pca = prcomp(t(imputedData))
  summary(pca)
  variance = round((pca$sdev^2 / sum(pca$sdev^2))*100, digits = 2)
  PCA_out = as.data.frame(pca$x) %>% 
    add_column(variance) %>% 
    rownames_to_column(var = "factorName")
  
  return(PCA_out)
}

# a) plot PCA with imputed observations
plotPCA_factors = factorPCA(imp) %>% 
  inner_join(meanResTime)
p = ggplot(plotPCA_factors, aes(x = PC1, y = PC2, label = factorName, fill = meanResTime, color = meanResTime))
p + geom_point(size=4, shape = 21, stroke = 0.6, color = "grey50") + 
  geom_text(nudge_y = -60, color = "grey40", size = 2) +
  theme_bw() +
  xlab(paste('PC1 (', plotPCA_factors $variance[1],'%)', sep = ""))+
  ylab(paste('PC2 (', plotPCA_factors $variance[2],'%)', sep = "")) +
  coord_fixed() +
  xlim(-270, 500) +
  ylim(-410, 150) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  guides(fill=guide_colorbar(title=paste("mean residence", "time [min]", sep = "\n"))) 
#ggsave("figures/panels/figS4/PCA_factors.pdf", width = 10, height = 7, units = "cm")
