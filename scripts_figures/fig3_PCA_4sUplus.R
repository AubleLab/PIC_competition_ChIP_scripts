rm(list = ls())

library(tidyverse)
library(DESeq2)

# upload count table, with filtered out low count genes
countTableFilt = read.csv("data/countCerevisiae_genes_minusStranded__lessThan12zeros.csv") %>% 
  column_to_rownames(var = "X")

# load metadata - keep only samples with added 4sU
metadata = read.csv("data/metadata.csv") %>% 
  mutate(normFactor = mappedReads_pombe/2000000) %>% 
  mutate(combo = paste0(fraction, "_", added_4sU, "4sU")) %>% 
  mutate(combo = as.factor(combo)) %>% 
  mutate(timeInGal = as.factor(timeInGal)) %>% 
  dplyr::filter(added_4sU == "+")
rownames(metadata) = metadata$sampleNames

countTableFilt = countTableFilt[, rownames(metadata)]

#PCA plot function
printPCA = function(x, ntop = 500, metadata){
  # the function prepares the data, computes the selected principal components and returns also variablity
  rv = rowVars(assay(x)) # get variance in each row
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca = prcomp(t(assay(x)[select,])) # PCA analysis
  
  #get the variance explained by ech of the components
  variance<-round((pca$sdev^2 / sum(pca$sdev^2))*100, digits = 2)
  
  PCA_out = as.data.frame(pca$x) %>% 
    add_column(variance) %>% 
    rownames_to_column(var = "sampleNames") %>% 
    left_join(metadata, by = "sampleNames")
  
  #summary(pca)
  return(PCA_out)
}

# get normalization factors
norm_matrix = matrix(rep(metadata$normFactor, nrow(countTableFilt)), byrow = T, nrow = nrow(countTableFilt))
colnames(norm_matrix) = colnames(countTableFilt)

# create DESeq2 object (to rlog normalize pombe-normalized counts)
dds = DESeqDataSetFromMatrix(countData=countTableFilt,colData=metadata, design=~timeInGal)

# add pombe normalization factors
normalizationFactors(dds) = norm_matrix


# transform the data
# rld = rlogTransformation(dds,blind=TRUE)
# save(rld, file = "data/rld_PCA_pombeNormalized_4sUonly.Rdata")
load("data/rld_PCA_pombeNormalized_4sUonly.Rdata")

# PCA
plotData = printPCA(rld, ntop = Inf, metadata)

# plot
colorPalette = c("#E29930", "#217CA3", "#32384D")
p = ggplot(plotData, aes(x=PC1, y=PC2, fill = fraction, color = timeInGal))
p + geom_point(size=6, alpha = 0.7, shape = 21, stroke = 1) +
  xlab(paste('PC1 (', plotData$variance[1],'%)', sep = "")) +
  ylab(paste('PC2 (', plotData$variance[2],'%)', sep = "")) +
  theme_bw() + 
  coord_fixed() + 
  scale_color_manual(values = c("white", "black"))+
  #xlim(-150, 250) +
  ylim(-30, 45) +
  scale_fill_manual(values = colorPalette)+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=12),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave("figures/panels/fig3/PCA.pdf", width = 13, height = 6, units = "cm")
