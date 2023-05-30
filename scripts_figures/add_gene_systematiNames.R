rm(list = ls())

#dd columns with systematic names

resTimes = read.csv("data/residence_times_all.csv") 

geneConvertTable = read.delim("/Users/lilcrusher/annotations/yeast/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)


resTimeEdit = resTimes %>% 
  left_join(geneConvertTable) %>% 
  distinct() %>% 
  dplyr::select(peakName, gene, geneSys, TBP:TFIIF)
write.csv(resTimeEdit, "data/residence_times_all.csv", quote = F, row.names = F)
