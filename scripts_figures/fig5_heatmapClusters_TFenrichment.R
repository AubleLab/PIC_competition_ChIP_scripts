rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ggrepel)

# upload heatmap clusters
allHeatmapGenes = read.csv("data/analysis/heatmap_clusters/heatmapClusters.csv")

# upload TF targets: assigned TFs to genes (only for regions withing -250:+100 bp from TSS)
load("data/YeastEpi_geneTargets.Rdata")


# do Fisher's test for each gene list vs all TF targets
for (clusterNum in levels(factor(allHeatmapGenes$cluster))){
  
  clusterGenes = allHeatmapGenes %>% 
    dplyr::filter(cluster == clusterNum)
  clusterGenes = clusterGenes$gene
  
  # go through the TF targets
  for (TF in names(geneTargets)){
    
    # get the targets for a give 
    TFgenes = geneTargets[[TF]]
    
    yes_yes = sum(clusterGenes %in% TFgenes)
    clusterYes_TFtargetNo = length(clusterGenes) - yes_yes
    TFtargetYes_clusterNo = length(TFgenes) - yes_yes
    
    leftover = length(unique(c(allHeatmapGenes$gene, TFgenes))) - yes_yes - clusterYes_TFtargetNo - TFtargetYes_clusterNo
    
    
    # create contingency table and o Fisher's exact test
    tbl = matrix(data = c(yes_yes, clusterYes_TFtargetNo, TFtargetYes_clusterNo, leftover), nrow = 2, byrow = TRUE)
    fisher = fisher.test(tbl, alternative = "greater")
    
    if (TF == names(geneTargets)[1] & clusterNum == "1"){
      resultTable = data.frame(cluster = clusterNum,
                               TF = TF,overlap = yes_yes, 
                               TFtargetYes_clusterNo = TFtargetYes_clusterNo, 
                               clusterYes_TFtargetNo = clusterYes_TFtargetNo, 
                               neither =  leftover, 
                               p = fisher$p.value)
    } else {
      helpTable = data.frame(cluster = clusterNum,
                                       TF = TF,overlap = yes_yes, 
                                       TFtargetYes_clusterNo = TFtargetYes_clusterNo, 
                                       clusterYes_TFtargetNo = clusterYes_TFtargetNo, 
                                       neither =  leftover, 
                                       p = fisher$p.value)
      
      resultTable = rbind(resultTable, helpTable)
    }
  }
  
}

# attach FDR
FDR = p.adjust(resultTable$p, method = "fdr")

resultTable$FDR = FDR 
resultTable$logFDR = -log10(resultTable$FDR)

# write down the results
write.csv(resultTable,"data/analysis/heatmap_clusters/TF_enrichment/allTFenrichment.csv")

# filte out significant and write down results
sigRes = resultTable %>% dplyr::filter(FDR<0.05) 
write.csv(sigRes,"data/analysis/heatmap_clusters/TF_enrichment/sigTFenrichment.csv")

# filter out significant and write down results
sigPRes = resultTable %>% dplyr::filter(p<0.05) 
write.csv(sigPRes,"data/analysis/heatmap_clusters/TF_enrichment/sigP_TFenrichment.csv")


