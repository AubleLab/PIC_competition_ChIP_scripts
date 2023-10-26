rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ggrepel)


# upload table with residence times and replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(c(gene, TFIIF)) %>% 
  drop_na() %>% 
  distinct()

plotResTime = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(TFIIF = ifelse(TFIIF == "<1", randNum, as.numeric(TFIIF))) %>% 
  dplyr::select(gene, TFIIF)



# upload residence time table - select TFIIF 
# and separate those into short-lived (5 min thresh.) or not
# assign short-lived labels
shortTFIIF = plotResTime %>% 
  mutate(shortTFIIF = ifelse(TFIIF <5, "yes", "no")) %>% 
  mutate(shortTFIIF = factor(shortTFIIF, levels = c("yes", "no")))

# upload TF targets
load("/Users/lilcrusher/annotations/yeast_DBFs/geneTargets.Rdata")


# do Fisher's test for each TF 
for (TF in names(geneTargets)){
  TFgenes = geneTargets[[TF]]
  
  # find out if a list is in the list of TF target genes
  testTable = shortTFIIF %>% 
    mutate(TFtarget = ifelse(gene %in% TFgenes, "yes", "no"))%>% 
    mutate(TFtarget = factor(TFtarget, levels = c("yes", "no")))
  
  # create contingency table and o Fisher's exact test
  tbl = table(testTable$TFtarget, testTable$shortTFIIF)
  fisher = fisher.test(tbl, alternative = "greater")
  
  if (TF == names(geneTargets)[1]){
    resultTable = data.frame(TF = TF,overlap = tbl[1,1], 
                             TFtargetYes_longTFFIFno = tbl[1,2], 
                             longTFFIFyes_TFtargetNo = tbl[2,1], 
                             neither = tbl[2,2], 
                             p = fisher$p.value)
  } else {
    helpTable = data.frame(TF = TF,overlap = tbl[1,1], 
                           TFtargetYes_longTFFIFno = tbl[1,2], 
                           longTFFIFyes_TFtargetNo = tbl[2,1], 
                           neither = tbl[2,2], 
                           p = fisher$p.value)
    
    resultTable = rbind(resultTable, helpTable)
  }
}

# attach FDR
FDR = p.adjust(resultTable$p, method = "fdr")

resultTable$FDR = FDR 
resultTable$logFDR = -log10(resultTable$FDR)
# remove PIC associated TFs from the results and plot
p = ggplot(resultTable %>% dplyr::filter(FDR<0.05) %>% 
             dplyr::filter(!str_detect(TF, "^Taf")) %>% 
             dplyr::filter(!str_detect(TF, "^CTD")) %>% 
             dplyr::filter(!str_detect(TF, "^Tf")) %>% 
             dplyr::filter(!str_detect(TF, "^Ccl1")) %>% 
             dplyr::filter(!str_detect(TF, "^Spt15")) %>% 
             dplyr::filter(!str_detect(TF, "^Ssl1")) %>% 
             dplyr::filter(!str_detect(TF, "^TAF")) %>% 
             dplyr::filter(!str_detect(TF, "^TOA1")) %>% 
             dplyr::filter(!str_detect(TF, "^TFA1")) %>% 
             dplyr::filter(!str_detect(TF, "^Kin28")) %>% 
             dplyr::filter(!str_detect(TF, "^Rpb")) %>% 
             dplyr::filter(!str_detect(TF, "^Bdf2")), 
           aes(y = reorder(TF, -FDR),x = logFDR, fill = logFDR))
p + geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9)) +
  xlab(expression(paste(-log[10],"(FDR)")))+
  ylab("")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sigRes = resultTable %>% dplyr::filter(FDR<0.05) #%>% 
  # dplyr::filter(!str_detect(TF, "^Taf")) %>% 
  # dplyr::filter(!str_detect(TF, "^CTD")) %>% 
  # dplyr::filter(!str_detect(TF, "^Tf")) %>% 
  # dplyr::filter(!str_detect(TF, "^Ccl1")) %>% 
  # dplyr::filter(!str_detect(TF, "^Spt15")) %>% 
  # dplyr::filter(!str_detect(TF, "^Ssl1")) %>% 
  # dplyr::filter(!str_detect(TF, "^TAF")) %>% 
  # dplyr::filter(!str_detect(TF, "^TOA1")) %>% 
  # dplyr::filter(!str_detect(TF, "^TFA1")) %>% 
  # dplyr::filter(!str_detect(TF, "^Kin28")) %>% 
  # dplyr::filter(!str_detect(TF, "^Rpb")) %>% 
  # dplyr::filter(!str_detect(TF, "^Bdf2"))

