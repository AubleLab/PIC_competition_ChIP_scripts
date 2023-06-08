rm(list = ls())

library(tidyverse)
library(annotateYeast)

# script for extracting additional fast sites

# first files with Hill fits
fitFiles = list.files("data/Hill_fits/", full.names = T, pattern = "*fullTableWithFits.csv")
for (i in fitFiles){
  fitTable = read.csv(i)
  
  if (i == fitFiles[1]){
    fullTableWithFits = fitTable
  } else (
    fullTableWithFits = rbind(fullTableWithFits, fitTable)
  )
}
#-------------------------------------------------------------------
# take the whole count table -> extract peaks -> 
# assign all peaks to the nearest genes ->
# filter peaks that are only within -250:100 distances from TSS
# remove t-RNAs
# if one gene has multiple peaks - keep the closer one
peaks = fullTableWithFits %>% 
  dplyr::select(peakName) %>%
  separate(peakName, c("chr", "start", "end"), sep="_") %>%
  mutate(start = as.integer(start)) %>%
  mutate(end = as.integer(end)) %>% 
  distinct()

peak_factor = fullTableWithFits %>% 
  dplyr::select(peakName, factorName)

# matching systematic names from www.yeastgenome.org:
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)

# actual peak filtering
distancesGenes = suppressWarnings(calcFeatureDist_aY(peaks, yeastGenes)) %>% 
  left_join(geneConvertTable) %>% 
  distinct() %>% 
  mutate(TSS = ifelse(strand == "+", startGene, endGene)) %>% 
  mutate(startPeak = ifelse(strand == "+", start, end)) %>% 
  mutate(endPeak = ifelse(strand == "+", end, start)) %>% 
  mutate(startDist = start - TSS) %>% 
  mutate(endDist = end - TSS) %>% 
  mutate(TSSdist = ifelse(startDist <= 0 & endDist >=0, 0,
                          ifelse(startDist < 0 & endDist <0, endDist, startDist))) %>% 
  mutate(absTSSdist = abs(TSSdist)) %>% 
  dplyr::select(-c(TSS, startPeak, endPeak, startDist, endDist)) %>% 
  filter(TSSdist > -250 & TSSdist < 100)%>% 
  filter(!str_detect(gene,"^t")) %>% 
  unite(peakName, chr:end, sep = "_") %>% 
  left_join(peak_factor) %>% 
  group_by(gene, factorName) %>% 
  dplyr::slice(which.min(absTSSdist)) %>% 
  ungroup() %>%
  dplyr::select(peakName, factorName) %>% 
  distinct() 

# rename genes back to the systematic names!
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6) %>% 
  dplyr::rename(geneName = gene) %>% 
  dplyr::rename(gene = geneSys)


myFits = fullTableWithFits %>% 
  inner_join(distancesGenes)%>% 
  left_join(geneConvertTable) %>% 
  mutate(gene = ifelse(!is.na(geneName), geneName, gene)) %>% 
  dplyr::select(-geneName)

# then upload reliable fits from the first estimate iteration - 
# these sites were reliably fitted and don't need to be fitted again
factorFiles = list.files("data/time_estimates_first_iter", 
                         pattern = ".csv", full.names = T)

for (i in factorFiles){
  factorTable = read.csv(i)
  
  splitName = unlist(strsplit(i, split = "_"))
  # find out if fast file or with time and find pit factor name
  fastOrTime = splitName[length(splitName)]
  factorName = toupper(splitName[(length(splitName) - 1)])
  
  # if fast replace time columns with random number below 1
  # in both add factor name
  if (fastOrTime == "fast.csv"){
    factorTableWithName = factorTable %>% 
      mutate(time = runif(nrow(factorTable), min=0.5, max=0.8)) %>% 
      mutate(factorName = factorName)
  } else if(fastOrTime == "ttfast.csv"){
    factorTableWithName = factorTable %>% 
      mutate(time = runif(nrow(factorTable), min=0, max=0.1)) %>% 
      mutate(factorName = factorName)
  } else {
    factorTableWithName = factorTable %>% 
      mutate(factorName = factorName)
  }
  
  # rbind the tables
  if (i == factorFiles[1]){
    finalFactorTable = factorTableWithName
  } else{
    finalFactorTable = rbind(finalFactorTable, factorTableWithName)
  }
  
}

# remove t-RNAs and make the table compatible with the "myFits" table
iterOneFits = finalFactorTable%>% 
  filter(!str_detect(gene,"^t'")) %>% 
  unite(peakName, c("chr", "start", "end"), sep = "_") %>% 
  dplyr::select(peakName, time) %>% 
  dplyr::rename(resTime = time)

rm(list = setdiff(ls(), c("myFits", "iterOneFits")))


# filter my fits for delta Km from nls model < 2
# join table with the fits from first iteration - keep only peaks that
# were not fitted in the first iteration
fastFits = myFits %>% 
  dplyr::select(peakName, gene, factorName, Km_western, Km_nls) %>% 
  mutate(delta_Km = Km_nls - Km_western) %>% 
  dplyr::filter(delta_Km < 2) %>% 
  full_join(iterOneFits) %>% 
  dplyr::filter(is.na(resTime)) %>% 
  dplyr::select(peakName, gene, factorName)
#write.csv(fastFits,"data/time_estimates_second_iter/additionalFastSites.csv", quote = F, row.names = F)




