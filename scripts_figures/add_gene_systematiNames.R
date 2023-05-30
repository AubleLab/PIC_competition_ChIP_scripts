rm(list = ls())

#dd columns with systematic gene names

resTimes = read.csv("data/residence_times_all.csv") 

# table from From: http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab :
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)


resTimeEdit = resTimes %>% 
  left_join(geneConvertTable) %>% 
  distinct() %>% 
  dplyr::select(peakName, gene, geneSys, TBP:TFIIF)
write.csv(resTimeEdit, "data/residence_times_all_withSysNames.csv", quote = F, row.names = F)
