rm(list = ls())
library(tidyverse)
#devtools::install_github("AubleLab/annotateYeast")
library(annotateYeast)
library(modelr)

# matching systematic names from www.yeastgenome.org:
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)

# ------ look at raw counts -------

# upload count tables- normalize to fall between 0-1
countFiles = list.files("data/normalized_count_tables", pattern = ".txt", full.names = T)

for (i in countFiles){
  countTableOrig = read.delim(i, sep = " ") 
  
  # get the last time point, so we can divide all the values by this column 
  # (this way) all values end at 1
  columnNames = colnames(countTableOrig)
  lastTimePoint = countTableOrig[,str_ends(columnNames, "120")]
  
  countTable = countTableOrig %>% 
    rownames_to_column("peakName") %>% 
    mutate(time120 = (lastTimePoint)) %>% 
    gather(-c(peakName, time120) ,key = "timeLabel", value = "normCount") %>% 
    mutate(normCount_ends1 = normCount / time120) %>% 
    group_by(peakName) %>% 
    mutate(normCount_01 = (normCount - min(normCount)) / (max(normCount) - min(normCount))) %>% 
    dplyr::select(-time120) %>% 
    mutate(time = readr::parse_number(timeLabel)) %>% 
    group_by(peakName) 
  
  splitName = unlist(strsplit(basename(i), split = "_"))
  # find out if fast file or with time and find pit factor name
  factorName = toupper(splitName[1])
  
  # if fast replace time columns with random number below 1
  # in both add factor name
  countTableWithName = countTable %>% 
    mutate(factorName = factorName)
  
  # rbind the tables
  if (i == countFiles[1]){
    finalCountTable = countTableWithName
  } else{
    finalCountTable = rbind(finalCountTable, countTableWithName)
  }
}
finalCountTable$factorName = factor(finalCountTable$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                           "TFIIF", "TFIIE"))

# upload western blots
finalWesternTable = read.csv("data/westerns.csv")


finalWesternTable = finalWesternTable %>% 
  dplyr::filter(time!=0)
finalWesternTable$factorName = factor(finalWesternTable$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                               "TFIIF", "TFIIE", "FKH1"))
# remove unnnecessary variables
rm(list = setdiff(ls(), c("geneConvertTable", "finalCountTable", "finalWesternTable")))


# ----- now try to fit all data - one factor at a time
# do fitting with the drc package - two parameter model on data that end at 1 / or 4-parameter fit
# with no restriction
# compare that to the fits incorporation the induction curve (Hill's coefficient) - either with
# saturation set to 1 or without saturation set

# upload fits for Westerns, from which we will get Hill coefficiennts
fitParamW = read.csv("data/westerns_fitParam.csv")

for (fitFactor in levels(finalCountTable$factorName)){
# updateLevels = levels(finalCountTable$factorName)[!levels(finalCountTable$factorName) %in% c("TBP", "TFIIA", "TFIIB")]
# for (fitFactor in updateLevels){  
  # filter a given factor and  get peaks
  factor_peaks = finalCountTable %>% 
    filter(factorName== fitFactor) %>% 
    dplyr::select(peakName) %>%
    separate(peakName, c("chr", "start", "end"), sep="_") %>%
    mutate(start = as.integer(start)) %>%
    mutate(end = as.integer(end)) %>% 
    distinct()
  
  # take the whole count table - assign all peaks to the nearest genes
  distancesGenes = suppressWarnings(calcFeatureDist_aY(factor_peaks, yeastGenes)) %>% 
    left_join(geneConvertTable) %>% 
    dplyr::select(gene, geneSys, chr, start, end) %>% 
    unite(peakName, c("chr", "start", "end")) %>% 
    distinct() %>% 
    dplyr::rename(geneName = gene)
  
  # combine residence times with gene names + add western blots
  factor_counts_all_wW = finalCountTable %>% 
    filter(factorName == fitFactor) %>% 
    left_join(distancesGenes) %>% 
    mutate(gene = ifelse(geneSys == "", geneName, geneSys)) %>% 
    dplyr::select(-c(geneName, geneSys)) %>% 
    left_join(finalWesternTable, by = c("factorName", "time")) %>% 
    distinct() %>% 
    mutate(geneLabel = factor(paste(peakName, gene, sep = ":")))
  
  # ---- do fitting with NLS, where n is set to n from western fits: ----
  # ---- but Vmax is also fitted - no to hard-set saturation point
  nlc = nls.control(maxiter = 1000, warnOnly = T)
  n = fitParamW %>% 
    dplyr::filter(factorName == fitFactor) %>% 
    dplyr::select(n_western)
  n = as.numeric(n)
  for (i in levels(factor_counts_all_wW$geneLabel)){
    
    tryCatch({
      nlsFit = nls(formula = normCount_ends1 ~ Vmax*time^(n) / (Km^(n) + time^(n)), 
                   factor_counts_all_wW %>% filter(geneLabel == i),
                   start=list(Km=40, Vmax = 1), 
                   control=nlc)
      coefs = coef(nlsFit)
      R2_nls = modelr::rsquare(nlsFit, factor_counts_all_wW %>% filter(geneLabel == i))
      
      tableNLS = data.frame(Km_nls = coefs["Km"], Vmax_nls = coefs["Vmax"], R2_nls, geneLabel = i)
      if (i == levels(factor_counts_all_wW$geneLabel)[1]){
        fitTableNLS_n = tableNLS
      } else {
        fitTableNLS_n = rbind(fitTableNLS_n, tableNLS)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
  }
  
  
  
  # add Km and n to the count table
  fullTableWithFits = factor_counts_all_wW %>% 
    ungroup() %>% 
    dplyr::select(peakName, factorName, gene, geneLabel) %>% 
    distinct() %>% 
    inner_join(fitTableNLS_n)
  
  # save the fits:
  #write.csv(fullTableWithFits, paste0("data/Hill_fits/",fitFactor,"_fullTableWithFits.csv"), row.names = F)
  
}


