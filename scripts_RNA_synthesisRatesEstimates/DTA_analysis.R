rm(list = ls())

library(DTA)
library(tidyverse)
# ====================================================================================================
# test 
data(Miller2011)
data(Miller2011dynamic)
# ---- DTA data
# Totals = Sc.datamat[,which(Sc.phenomat[,"fraction"]=="T")]
# Total = apply(log(Totals),1,median)
# plotsfkt = function(){
#   par(mar = c(5,4,4,2)+0.1+1)
#   par(mai = c(1.1,1.1,1.3,0.7))
#   hist(Total,breaks = seq(0,ceiling(max(Total)),1/4),
#          cex.main=1.5,cex.lab=1.25,cex.axis=1.125,
#          main="Histogram of log(Total)",
#          xlab="gene-wise median of total samples")
#   hist(Total[Total >= 5],breaks = seq(0,ceiling(max(Total)),1/4),
#          col = "#08306B",add = TRUE)
#   }
# DTA.plot.it(filename = "figures/DTA_test/histogram_cut_off",plotsfkt = plotsfkt,
#               saveit = TRUE,notinR = TRUE)
# 
# dev.off()
# 
# res = DTA.estimate(Sc.phenomat,Sc.datamat,Sc.tnumber,
#                    ccl = 150,mRNAs = 60000,reliable = Sc.reliable,
#                    condition = "real_data",save.plots = TRUE,notinR = TRUE,folder = "figures/DTA_test")
# dev.off()
# 
# 
# # --- let's have a look at their data ---
# 
# plotTheirCounts = as.data.frame(Sc.datamat) %>%
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(as.data.frame(Sc.phenomat)) %>% 
#   mutate(logIntensity = log(normCount)) %>% 
#   mutate(fraction = factor(fraction, levels = c("L", "U", "T")))
# 
# p = ggplot(plotTheirCounts, aes(logIntensity, fill = fraction))
# p + geom_histogram(bins = 50) +
#   theme_bw() +
#   facet_wrap(~name)
# ggsave("figures/DTA_test/my_plots/intensity_histogram.pdf", width = 8, height = 5)
#   
# # get mean value of the log transformed intensities
# plotThierSummarized = plotTheirCounts %>% 
#   group_by(gene, fraction) %>% 
#   summarise(meanLog = mean(logIntensity)) 
# 
# p = ggplot(plotThierSummarized, aes(meanLog, fill = fraction))
# p + geom_histogram(bins = 50) +
#   theme_bw() +
#   facet_wrap(~fraction, ncol = 1)
# ggsave("figures/DTA_test/my_plots/intensity_histogram_meanVal.pdf", width = 3, height = 5)

# ====================================================================================================
# ========= our data =======
# ====================================================================================================
countTable = read.csv("countTables/filtered_countTables/countCerevisiae_genes_minusStranded__ThiouracilYESonly_pombeNormalized.csv") %>% 
  column_to_rownames(var = "X")

# load metadata and filter out the samples where thiouracil was not added
metadata = read.csv("metadata.csv") %>% 
  filter(thiouracil == "yes") %>% 
  dplyr::rename(name = sampleNames) %>% 
  mutate(time = 6) %>% 
  mutate(fraction = ifelse(sampleType == "biotin", "L", 
                           ifelse(sampleType == "total", "T", "U")))
rownames(metadata) = metadata$name
metadata$nr = rep(c(1,2,3,4), 3)

# filter count non-thiouracil samples from countTable and convert to matrix
countTable = countTable[, rownames(metadata)]

#let's have a look at the normalized counts across each condition
plotCounts = countTable %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
  left_join(metadata) %>% 
  mutate(logCounts = log(normCount))

p = ggplot(plotCounts, aes(logCounts, fill = sampleType))
p + facet_wrap(~name) +
  theme_bw() +
  geom_histogram(bins = 50)
ggsave("figures/DTA/normCount_histogram.pdf", width = 8, height = 5)

plotCountsSummarised = plotCounts %>%
  group_by(gene, sampleType) %>% 
  summarise(meanLog = mean(logCounts)) 
p = ggplot(plotCountsSummarised, aes(meanLog, fill = sampleType))
p + facet_wrap(~sampleType, ncol = 1) +
  theme_bw() +
  geom_histogram(bins = 50)
ggsave("figures/DTA/normCount_histogram_meanVal.pdf", width = 3.5, height = 5)

# it seems like there cannot be any zeros in the matrix
zeroSum = rowSums(countTable == 0)
zeroSum = zeroSum[zeroSum == 0]

countTable_nonZero = countTable[names(zeroSum),]


# Sc.datamat = countTable_DTA
# Sc.phenomat = metadata_DTA
# Sc.tnumber - from the package
# Sc. reliable = reliable_DTA
reliable_DTA = rownames(countTable_nonZero)
metadata_DTA = as.matrix(metadata)
countTable_DTA = as.matrix(countTable_nonZero)

# res_DTA = DTA.estimate(phenomat = metadata_DTA,
#                        datamat = as.matrix(countTable),
#                        tnumber = Sc.tnumber,
#                        check = T,
#                        ccl = 150,
#                        mRNAs = 60000,
#                        reliable = reliable_DTA, 
#                        LtoTratio = 0.1,
#                        condition = "real_data",
#                        save.plots = TRUE,notinR = TRUE,folder = "figures/DTA")

# use method bias to estimate LtoT ratio since default failed
res_DTA = DTA.estimate(phenomat = metadata_DTA,
                       datamat = countTable_DTA,
                       tnumber = Sc.tnumber, 
                       check = T,
                       ccl = 150,
                       mRNAs = 60000,
                       reliable = reliable_DTA,
                       condition = "real_data",save.plots = TRUE,notinR = TRUE,folder = "figures/DTA", 
                       ratiomethod = "bias")

resList = res_DTA$`6`

# extract results
decay = resList$dr
synthesis = resList$sr
halfLife = resList$hl
genes = names(synthesis)

resTable = data.frame(gene = genes, decay, synthesis, halfLife)
# find out how many genes have reliable estimated
resTableFilt = resTable %>% filter(!is.nan(synthesis))
write.csv(resTable, paste0("/Users/lilcrusher/yeast_RNA/DTA_results/results/res_biasMethod.csv"))


# lets' try setting different LtoTratios
for (i in seq(0.01, 0.1, by = 0.01)){
  
  iName = str_replace(as.character(i), pattern = "[.]", replacement = "_")
  
  figDir = paste0("/Users/lilcrusher/yeast_RNA/DTA_results/figures/LtoT_",i)
  
  ifelse(file.exists(figDir), NA, dir.create(figDir))
  
  res_DTA = DTA.estimate(phenomat = metadata_DTA,
                         datamat = countTable_DTA,
                         tnumber = Sc.tnumber, 
                         check = T,
                         ccl = 150,mRNAs = 60000,reliable = reliable_DTA,
                         condition = "real_data",save.plots = TRUE,notinR = TRUE,folder = figDir, 
                         LtoTratio = i)
  
  resList = res_DTA$`6`
  
  # extract results
  decay = resList$dr
  synthesis = resList$sr
  halfLife = resList$hl
  genes = names(synthesis)
  
  
  resTable = data.frame(gene = genes, decay, synthesis, halfLife)
  
  
  write.csv(resTable, paste0("/Users/lilcrusher/yeast_RNA/DTA_results/results/res_LtoT_", iName, ".csv"))
}

# # ====================================================================================================
# # ===== add +1 to the orig count table to eliminate 0 counts =====
# # ====================================================================================================
# countTable_plusOne = countTable + 1
# 
# #let's have a look at the normalized counts across each condition
# plotCounts_plusOne = countTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   mutate(logCounts = log(normCount))
# 
# p = ggplot(plotCounts_plusOne, aes(logCounts, fill = sampleType))
# p + facet_wrap(~name) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA/normCount_add1_histogram.pdf", width = 8, height = 5)
# 
# plotCounts_plusOneSummarised = plotCounts_plusOne %>%
#   group_by(gene, sampleType) %>% 
#   summarise(meanLog = mean(logCounts)) 
# p = ggplot(plotCounts_plusOneSummarised, aes(meanLog, fill = sampleType))
# p + facet_wrap(~sampleType, ncol = 1) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA/normCount_add1_histogram_meanVal.pdf", width = 3.5, height = 5)
# 
# 
# # ------- remove counts that have mean <2 > 1000 at least in one group
# meanVals = countTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   group_by(gene, sampleType) %>% 
#   summarise(meanVal = mean(normCount)) %>% 
#   pivot_wider(names_from = "sampleType", values_from = "meanVal") %>% 
#   column_to_rownames(var = "gene")
# 
# lessOne = rowSums(meanVals < 2 | meanVals > 1000)
# filtLessOne = lessOne[lessOne == 0]
# 
# # filtered count table
# filtTable_plusOne = countTable_plusOne[names(filtLessOne),]
# 
# plotCountsFilt_plusOne = filtTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   mutate(logCounts = log(normCount))
# 
# p = ggplot(plotCountsFilt_plusOne, aes(logCounts, fill = sampleType))
# p + facet_wrap(~name) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA/normCount_histogram__filtered.pdf", width = 8, height = 5)
# 
# plotCountsSummarisedFilt = plotCountsFilt_plusOne %>%
#   group_by(gene, sampleType) %>% 
#   summarise(meanLog = mean(logCounts)) 
# p = ggplot(plotCountsSummarisedFilt, aes(meanLog, fill = sampleType))
# p + facet_wrap(~sampleType, ncol = 1) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# #ggsave("figures/DTA/normCount_histogram_meanVal__filtered.pdf", width = 3.5, height = 5)
# 
# 
# reliable_DTA_plusOne = rownames(filtTable_plusOne)
# countTable_DTA_plusOne = as.matrix(filtTable_plusOne)
# 
# res_DTA = DTA.estimate(phenomat = metadata_DTA,
#                        datamat = countTable_DTA_plusOne,
#                        tnumber = Sc.tnumber, 
#                        check = T,
#                        ccl = 150,mRNAs = 60000,
#                        reliable = reliable_DTA_plusOne,
#                        condition = "real_data",save.plots = TRUE,
#                        notinR = TRUE,folder = "DTA_results/figures/counts_plusOne")
# 
# resList = res_DTA$`6`
# 
# # extract results
# decay = resList$dr
# synthesis = resList$sr
# halfLife = resList$hl
# genes = names(synthesis)
# 
# resTable = data.frame(gene = genes, decay, synthesis, halfLife)
# 
# write.csv(resTable, "/Users/lilcrusher/yeast_RNA/DTA_results/results/counts_plusOne.csv")

# ====================================================================================================
# ======= let's now compare different approaches =========
# ====================================================================================================
rm(list = ls())
library(GGally)

resPath = "/Users/lilcrusher/yeast_RNA/DTA_results/results"
resFiles = list.files(resPath, full.names = T)


for (i in resFiles){
  
  fileName = unlist(strsplit(basename(i), split = ".", fixed = T))[1]
  
  helpRes = read.csv(i) %>% 
    mutate(fileName = fileName)
  
  if (i == resFiles[1]){
    resTable = helpRes
  } else {
    resTable = rbind(resTable, helpRes)
  }
  
}

# compare synthesis rates with different settings
synthesisRates = resTable %>% 
  select(gene, synthesis, fileName) %>% 
  pivot_wider(values_from = synthesis, names_from = fileName)

ggpairs(synthesisRates, columns = 2: 12) +
  theme_bw() + xlim(0, 2000) + ylim(0, 2000)
ggsave("DTA_results/compare_synthesis_cropped_2000.png", width = 15, heigh = 15)


# ====================================================================================================
# ====== now let's compare our synthesis rates LtoT=0.09 to the published data =================
# ====================================================================================================
rm(list = ls())

DTAres = read.csv("/Users/lilcrusher/yeast_RNA/DTA_results/results/res_LtoT_0_09.csv") %>% 
  select(gene, synthesis) %>% 
  mutate(synthesisPerMin = synthesis/150)
TR = read.delim("/Users/lilcrusher/competitionChIP/transcription_rates/transcription_rates_galactose.txt") %>% 
  select(-gene_name) %>% 
  dplyr::rename(gene = ORF) %>% 
  inner_join(DTAres) %>% 
  drop_na()

# crop outliers for plotting purposes
r = cor(TR$TR_galactose, TR$synthesisPerMin, method = "spearman")
p = ggplot(TR%>% 
             mutate(TR_galactose = ifelse(TR_galactose > 0.25, 0.25, TR_galactose)) %>% 
             mutate(synthesisPerMin = ifelse(synthesisPerMin > 2.5, 2.5, synthesisPerMin))
           , aes(x = synthesisPerMin, y = TR_galactose))
p + geom_point() +
  theme_bw() + 
  #ylim(0,0.25) +
  ggtitle(paste0("Spearman ro = ", round(r, digits = 2))) +
  xlab("our data (mRNA per cell per min)") +
  ylab("published data")
ggsave("DTA_results/compareSynthesisRates_LtoT009_vs_publishedData.pdf", width = 5, height = 5.5)

# ====================================================================================================
# ====== now let's compare our synthesis rates to the published data =================
# ====================================================================================================
rm(list = ls())

DTAres = read.csv("DTA_results/results/res_biasMethod.csv") %>% 
  select(gene, synthesis) %>% 
  mutate(synthesisPerMin = synthesis/150)
TR = read.delim("/Users/lilcrusher/competitionChIP/transcription_rates/transcription_rates_galactose.txt") %>% 
  select(-gene_name) %>% 
  dplyr::rename(gene = ORF) %>% 
  inner_join(DTAres) %>% 
  drop_na()

# crop outliers for plotting purposes
r = cor(TR$TR_galactose, TR$synthesisPerMin, method = "spearman")
p = ggplot(TR%>% 
             mutate(TR_galactose = ifelse(TR_galactose > 0.25, 0.25, TR_galactose)) %>% 
             mutate(synthesisPerMin = ifelse(synthesisPerMin > 1, 1, synthesisPerMin))
           , aes(x = synthesisPerMin, y = TR_galactose))
p + geom_point() +
  theme_bw() + 
  #ylim(0,0.25) +
  ggtitle(paste0("Spearman ro = ", round(r, digits = 2))) +
  xlab("our data (mRNA per cell per min)") +
  ylab("published data") 
ggsave("DTA_results/compareSynthesisRates_biasMethod_vs_publishedData.pdf", width = 5, height = 5.5)

###--------create a FINAL table with synthesis rated - mRNA per cel per minute------
rm(list = ls())
DTAres = read.csv("DTA_results/results/res_biasMethod.csv") %>% 
  select(gene, synthesis) %>% 
  mutate(synthesis_perCell_perMin = synthesis/150) %>% 
  filter(!is.na(synthesis))
write.csv(DTAres,"/Users/lilcrusher/yeast_RNA/DTA_results/FINAL_results_to_use/res_biasMethod.csv", 
          row.names = F)

p = ggplot(DTAres, aes(x = synthesis_perCell_perMin))
p + geom_histogram(bins = 100, color = "darkred", fill = "darkred", alpha = 0.5) +
   theme_bw() +
  xlim(0, 1)
ggsave("figures/synthesisRates_biasMethod_histogram.pdf", width = 4, height = 3)

# ====================================================================================================
# ====== now let's compare our synthesis rates to mean values from biotin samples =================
# ====================================================================================================
rm(list = ls())

DTAres = read.csv("DTA_results/FINAL_results_to_use/res_biasMethod.csv") %>% 
  select(gene, synthesis)
biotin = read.csv("/Users/lilcrusher/yeast_RNA/countTables/mean_counts_sampleTypes/biotin.csv") %>% 
  dplyr::rename(gene = ORF) %>% 
  inner_join(DTAres)%>% 
  drop_na()

r = cor(biotin$normMeanCount, biotin$synthesis, method = "spearman")
p = ggplot(biotin, aes(x = synthesis, y = normMeanCount))
p + geom_point() +
  theme_bw() + 
  ggtitle(paste0("Spearman ro = ", round(r, digits = 2))) +
  xlab("synthesis") +
  ylab("normalized mean count - biotin data") +
  xlim(0,500) +
  ylim(0, 10000)
ggsave("DTA_results/compareSynthesisRates_vs_meanNormCount_biotin.pdf", width = 5, height = 5.5)


### following FAILED with new count tables
# # ====================================================================================================
# # ====== let's try to get results from un normalized data =================
# # ====================================================================================================
# rm(list = ls())
# 
# data(Miller2011)
# 
# countTable = read.csv("countTables/countCerevisiae_genes.csv") %>% 
#   column_to_rownames(var = "X")
# 
# editColumns = unlist(strsplit(colnames(countTable), split = "_", fixed = T))[seq(1, 2*ncol(countTable), by = 2)]
# colnames(countTable) = editColumns
# 
# # load metadata and filter out the samples where thiouracil was not added
# metadata = read.csv("metadata.csv") %>% 
#   filter(thiouracil == "yes") %>% 
#   dplyr::rename(name = sampleNames) %>% 
#   mutate(time = 6) %>% 
#   mutate(fraction = ifelse(sampleType == "biotin", "L", 
#                            ifelse(sampleType == "total", "T", "U")))
# rownames(metadata) = metadata$name
# metadata$nr = rep(c(1,2,3,4), 3)
# 
# # filter count non-thiouracil samples from countTable and convert to matrix
# countTable = countTable[, rownames(metadata)]
# 
# #let's have a look at the normalized counts across each condition
# plotCounts = countTable %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   mutate(logCounts = log(normCount))
# 
# p = ggplot(plotCounts, aes(logCounts, fill = sampleType))
# p + facet_wrap(~name) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_histogram.pdf", width = 8, height = 5)
# 
# plotCountsSummarised = plotCounts %>%
#   group_by(gene, sampleType) %>% 
#   summarise(meanLog = mean(logCounts)) 
# p = ggplot(plotCountsSummarised, aes(meanLog, fill = sampleType))
# p + facet_wrap(~sampleType, ncol = 1) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_histogram_meanVal.pdf", width = 3.5, height = 5)
# 
# # ====================================================================================================
# # ===== add +1 to the orig count table to eliminate 0 counts =====
# countTable_plusOne = countTable + 1
# 
# #let's have a look at the normalized counts across each condition
# plotCounts_plusOne = countTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   mutate(logCounts = log(normCount))
# 
# p = ggplot(plotCounts_plusOne, aes(logCounts, fill = sampleType))
# p + facet_wrap(~name) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_add1_histogram.pdf", width = 8, height = 5)
# 
# plotCounts_plusOneSummarised = plotCounts_plusOne %>%
#   group_by(gene, sampleType) %>% 
#   summarise(meanLog = mean(logCounts)) 
# p = ggplot(plotCounts_plusOneSummarised, aes(meanLog, fill = sampleType))
# p + facet_wrap(~sampleType, ncol = 1) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_add1_histogram_meanVal.pdf", width = 3.5, height = 5)
# 
# 
# # ------- remove counts that have mean <2 > 1000 at least in one group
# meanVals = countTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   group_by(gene, sampleType) %>% 
#   summarise(meanVal = mean(normCount)) %>% 
#   pivot_wider(names_from = "sampleType", values_from = "meanVal") %>% 
#   column_to_rownames(var = "gene")
# 
# lessOne = rowSums(meanVals < 2 | meanVals > 1000)
# filtLessOne = lessOne[lessOne == 0]
# 
# # filtered count table
# filtTable_plusOne = countTable_plusOne[names(filtLessOne),]
# 
# plotCountsFilt_plusOne = filtTable_plusOne %>% 
#   rownames_to_column(var = "gene") %>% 
#   pivot_longer(!gene, values_to = "normCount", names_to = "name") %>% 
#   left_join(metadata) %>% 
#   mutate(logCounts = log(normCount))
# 
# p = ggplot(plotCountsFilt_plusOne, aes(logCounts, fill = sampleType))
# p + facet_wrap(~name) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_histogram__filtered.pdf", width = 8, height = 5)
# 
# plotCountsSummarisedFilt = plotCountsFilt_plusOne %>%
#   group_by(gene, sampleType) %>% 
#   summarise(meanLog = mean(logCounts)) 
# p = ggplot(plotCountsSummarisedFilt, aes(meanLog, fill = sampleType))
# p + facet_wrap(~sampleType, ncol = 1) +
#   theme_bw() +
#   geom_histogram(bins = 50)
# ggsave("figures/DTA_notNorm/normCount_histogram_meanVal__filtered.pdf", width = 3.5, height = 5)
# 
# metadata_DTA = as.matrix(metadata)
# reliable_DTA_plusOne = rownames(filtTable_plusOne)
# countTable_DTA_plusOne = as.matrix(filtTable_plusOne)
# 
# res_DTA = DTA.estimate(phenomat = metadata_DTA,
#                        datamat = countTable_DTA_plusOne,
#                        tnumber = Sc.tnumber, 
#                        LtoTratio = 0.01,
#                        check = T,
#                        ccl = 150,mRNAs = 60000,
#                        reliable = reliable_DTA_plusOne,
#                        condition = "real_data",save.plots = TRUE,
#                        notinR = TRUE,folder = "DTA_results_notNorm/figures/counts_plusOne")
# 
# resList = res_DTA$`6`
# 
# # extract results
# decay = resList$dr
# synthesis = resList$sr
# halfLife = resList$hl
# genes = names(synthesis)
# 
# resTable = data.frame(gene = genes, decay, synthesis, halfLife)
# 
# write.csv(resTable, "/Users/lilcrusher/yeast_RNA/DTA_results_notNorm/results/counts_plusOne.csv")
# 
# # ====================================================================================================
# # ======= let's now compare different approaches =========
# # ====================================================================================================
# rm(list = ls())
# 
# norm = read.csv("DTA_results/results/counts_plusOne.csv") %>% 
#   mutate(fileName = "norm")
# notNorm = read.csv("DTA_results_notNorm/results/counts_plusOne.csv")%>% 
#   mutate(fileName = "notNorm")
# resTable = rbind(norm, notNorm)
# 
# # compare synthesis rates with different settings
# synthesisRates = resTable %>% 
#   select(gene, synthesis, fileName) %>% 
#   pivot_wider(values_from = synthesis, names_from = fileName)
# 
# p = ggplot(synthesisRates, aes(x = norm, y = notNorm ))
# p +geom_point() +
#   theme_bw()
# ggsave("DTA_results_notNorm/compare_synthesis_norm_vs_notNorm.pdf", width = 5, heigh = 5)
# 
# # ====================================================================================================
# # ====== now let's compare these synthesis rates to the published data =================
# 
# rm(list = ls())
# 
#   DTAres = read.csv("DTA_results_notNorm/results/counts_plusOne.csv") %>% 
#   select(gene, synthesis)
# TR = read.delim("/Users/lilcrusher/competitionChIP/transcription_rates/transcription_rates_galactose.txt") %>% 
#   select(-gene_name) %>% 
#   dplyr::rename(gene = ORF) %>% 
#   inner_join(DTAres) %>% 
#   drop_na()
# 
# r = cor(TR$TR_galactose, TR$synthesis, method = "spearman")
# p = ggplot(TR%>% 
#              mutate(TR_galactose = ifelse(TR_galactose > 0.25, 0.25, TR_galactose)) %>% 
#              mutate(synthesis = ifelse(synthesis > 5, 5, synthesis))
#            , aes(x = synthesis, y = TR_galactose))
# p + geom_point() +
#   theme_bw() + 
#   ylim(0,0.25) +
#   ggtitle(paste0("Spearman ro = ", round(r, digits = 2))) +
#   xlab("our data") +
#   ylab("published data")
# ggsave("DTA_results_notNorm/compareSynthesisRates_vs_publishedData.pdf", width = 5, height = 5.5)
# 
# 
# 
# 
# 
# 
# 
# 
