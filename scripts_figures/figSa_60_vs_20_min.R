rm(list = ls())

library(tidyverse)
library(DESeq2)

# upload count table, with filtered out low count genes
countTableFilt = read.csv("data/countCerevisiae_genes_minusStranded__lessThan12zeros.csv") %>% 
  column_to_rownames(var = "X")

# load metadata - keep only samples with added 4sU
metadata = read.csv("data/metadata.csv") %>% 
  mutate(normFactor = mappedReads_pombe/2000000) %>% 
  mutate(fraction = as.factor(fraction)) %>% 
  mutate(timeInGal = as.factor(timeInGal)) %>% 
  dplyr::filter(added_4sU == "+")
rownames(metadata) = metadata$sampleNames

countTableFilt = countTableFilt[, rownames(metadata)]

# define function to get data to polt MA plot
makeMAplotTable = function(dds, IHW=F){
  #get normalized means and transform to log10
  if (IHW == T){
    res = as.data.frame(results(dds, filterFun=ihw))
  } else {
    res = as.data.frame(results(dds))
  }
  
  
  res$label = rep('non-significant', nrow(res))
  res$label[res$pvalue <= 0.05]="significant_p"
  res$label[res$padj <= 0.05]="significant_padj"
  
  
  res$label2 <- rep('non_significant', dim(res)[1])
  res$label2[res$padj <= 0.05]="significant_padj"
  
  res$log10baseMean = log10(res$baseMean)
  #get p-vales and transform to log10 scale
  res$logPval <- -10*log10(res$pvalue)
  
  return(res)
}

# perform DESeq2 analysis for every fraction with added thiouracil - compare
# 20 vs 60 min time points
for (i in levels(metadata$fraction)){
  
  subMetadata = metadata %>% 
    filter(fraction == i)
  
  subTable = countTableFilt[, rownames(subMetadata)]
  
  # get normalization factors
  norm_matrix = matrix(rep(subMetadata$normFactor, nrow(subTable)), byrow = T, nrow = nrow(subTable))
  colnames(norm_matrix) = colnames(subTable)
  
  # create DESeq2 object
  dds = DESeqDataSetFromMatrix(countData=subTable,colData=subMetadata, design=~timeInGal)
  
  # add pombe normalization factors
  normalizationFactors(dds) = norm_matrix
  dds$timeInGal <-relevel(dds$timeInGal, ref="20")
  dds = DESeq(dds)
  
  
  res = results(dds)
  res = as.data.frame(res)
  res = res[order(res$padj),]
  #write.csv(res,paste0("DESeq2/results/",i,".csv"), quote = F)
  
  MAplot = res %>% 
    mutate(label = ifelse(is.na(padj) | padj > 0.05, "non_sig", "sig")) %>% 
    mutate(fraction = i)
  
  sigUp = res %>% 
    filter(padj < 0.05 & log2FoldChange > 0) %>% 
    nrow()
  
  sigDown = res %>% 
    filter(padj < 0.05 & log2FoldChange < 0) %>% 
    nrow()
  upDown = data.frame(sigUp, sigDown) %>% 
    mutate(fraction = i)
  
  if (i == levels(metadata$fraction)[1]){
    plotTable = MAplot
    plotText = upDown
  } else{
    plotTable = rbind(plotTable, MAplot)
    plotText = rbind(plotText, upDown)
  }

}

plotTable$fraction = factor(plotTable$fraction, levels = c("L", "U", "T"))
plotText$fraction = factor(plotText$fraction, levels = c("L", "U", "T"))

# or for results with independent hypothesis testing padj:
#plotTable = makeMAplotTable(dds, IHW = T)
#make custom MA plot with red p-vals
p = ggplot(plotTable, aes(x = baseMean, y = log2FoldChange))
p + geom_point(color = "grey60", alpha = 0.3)+
  geom_point(data = plotTable %>% filter(label == "sig" & log2FoldChange < 0), color = "#E64A45")+ 
  geom_point(data = plotTable %>% filter(label == "sig" & log2FoldChange >0), color = "#96C0CE") +
  scale_x_continuous(trans = "log10", labels = function(x) format(x, scientific = F)) +
  theme_classic(base_size = 10)+
  #ylim(-0.7, 0.2) +
  geom_hline(yintercept = 0, lwd = 1) + 
  xlab("mean of normalized counts") + 
  ylab(expression(paste(log[2],"(fold-change)", " 60 vs.20 min"))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))  +
  geom_text(data = plotText, aes(x = 20000, y = 7, label = paste0("up = ", sigUp)), color = "#96C0CE") +
  geom_text(data = plotText, aes(x = 20000, y = -7, label = paste0("down = ", sigDown)), color = "#E64A45") +
  ylim(-7, 7) +
  facet_wrap(~fraction)
ggsave("figures/panels/figSa/MAplots_60_vs20_min.pdf", width = 15, height = 6, units = "cm")


