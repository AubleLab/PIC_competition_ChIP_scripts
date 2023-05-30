rm(list = ls())

library(tidyverse)
library(missMDA)

set.seed(42)

# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-peakName) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

geneTable = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime))) %>% 
  dplyr::select(gene, factorName, resTimesNum) %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  column_to_rownames("gene")


# impute the missing data with missMDA package
nb = estim_ncpPCA(geneTable,method.cv = "Kfold", verbose = FALSE)
res.comp = imputePCA(geneTable, ncp = nb$ncp)
imp = res.comp$completeObs

# save the imputed values
write.csv(imp, "data/imputed_missing_values/imputed_residence_times.csv")


#### keep only values that have residence time > 1
rm(list = ls())

library(tidyverse)
library(missMDA)

set.seed(42)

# upload table replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-peakName) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

geneTable = resTimes %>% 
  dplyr::filter(resTime != "<1") %>% 
  mutate(resTimesNum =  as.numeric(resTime)) %>% 
  dplyr::select(gene, factorName, resTimesNum) %>% 
  pivot_wider(names_from = factorName, values_from = resTimesNum) %>% 
  column_to_rownames("gene")


# impute the missing data with missMDA package
nb = estim_ncpPCA(geneTable,method.cv = "Kfold", verbose = FALSE)
res.comp = imputePCA(geneTable, ncp = nb$ncp)
imp = res.comp$completeObs

# save the imputed values
write.csv(imp, "data/imputed_missing_values/imputed_residence_times_fastSitesExcluded.csv")











