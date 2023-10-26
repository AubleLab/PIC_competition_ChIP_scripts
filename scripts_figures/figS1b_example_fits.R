rm(list = ls())
library(tidyverse)

# 1) =====================================================================================
# upload normalized count tables with HA/Myc ratios (I will take data in this raw form and normalize in a way that
# all end at 1)
countFiles = list.files("data//normalized_count_tables", pattern = ".txt", full.names = T)

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
    dplyr::select(-time120) %>% 
    mutate(time = readr::parse_number(timeLabel)) %>% 
    group_by(peakName) 
  
  splitName = unlist(strsplit(basename(i), split = "_"))
  # find out if fast file or with time and find fit factor name
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
rm(list = setdiff(ls(), "finalCountTable"))

# 2) =====================================================================================
# upload fits
# first upload my fit files
fitFiles = list.files("data/Hill_fits/", full.names = T, pattern = "*fullTableWithFits.csv")
for (i in fitFiles){
  fitTable = read.csv(i)
  
  if (i == fitFiles[1]){
    fullTableWithFits = fitTable
  } else (
    fullTableWithFits = rbind(fullTableWithFits, fitTable)
  )
}

fits = fullTableWithFits %>% 
  dplyr::select(peakName, Km_nls, Vmax_nls)

rm(list = setdiff(ls(), c("finalCountTable", "fits")))

# 3) =====================================================================================
# -- upload western blots (we use n from westerns)
western = read.csv("data/westerns.csv")
western$factorName = factor(western$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                           "TFIIF", "TFIIE"))

westernFits = read.csv("data/westerns_fitParam.csv")

# 4) =====================================================================================
# upload table with residence times replace <1 min with a random number between 0-1
# make class labels and select one random entry for each class
resTimes = read.csv("data/residence_times_all.csv") %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

countTimeGroups = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime)))%>% 
  mutate(timeGroup = ifelse(resTimesNum < 1, "<1 min",
                            ifelse(resTimesNum < 10, "1-10 min", ">10 min"))) 
# sampleTable = countTimeGroups %>% 
#   group_by(factorName, timeGroup) %>% slice_sample(n=1) %>% 
#   dplyr::select(peakName, gene, factorName, resTime, timeGroup) %>% 
#   inner_join(westernFits) %>% 
#   inner_join(fits)

# NOTE - load selected genes for plotting or use commented code above 
# that will take random numbers
examples = read.csv("data/genes_for_fits/examples.csv")
sampleTable = countTimeGroups %>% 
  inner_join(examples) %>% 
  dplyr::select(peakName, gene, factorName, resTime, timeGroup) %>% 
  inner_join(westernFits) %>% 
  inner_join(fits)

# reconstruct fits
x = seq(1, 125)
for (i in seq(1,nrow(sampleTable))){
  
  Vmax = sampleTable$Vmax_nls[i]
  Km = sampleTable$Km_nls[i]
  n = sampleTable$n_western[i]
  
  y =(Vmax*x^(n) / (Km^(n) + x^(n)))/Vmax
  
  fitTable = data.frame(time = x, y, peakName = rep(sampleTable$peakName[i], length(y))) %>% 
    left_join(sampleTable)
  if (i == 1){
    plotFitTable = fitTable
  } else {
    plotFitTable = rbind(plotFitTable, fitTable)
  }
}

# create a final table for plotting
plotSample = sampleTable %>% 
  inner_join(finalCountTable) %>% 
  dplyr::select(-c(normCount, timeLabel)) %>% 
  inner_join(western) %>% 
  inner_join(fits) %>% 
  full_join(plotFitTable)
plotSample$timeGroup = factor(plotSample$timeGroup, 
                              levels = c("<1 min", "1-10 min", ">10 min"))

plotSample$factorName = factor(plotSample$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                 "TFIIF", "TFIIE"))

textLabels = plotSample %>% 
  dplyr::select(gene, factorName, resTime) %>% 
  distinct() %>% 
  inner_join(countTimeGroups) %>% 
  mutate(resTime = ifelse(resTime == "<1", resTime, as.character(round(resTimesNum, digits = 1)))) %>% 
  mutate(myLabels = paste0(gene, ": ", resTime))
textLabels$timeGroup = factor(textLabels$timeGroup, 
                              levels = c("<1 min", "1-10 min", ">10 min"))

textLabels$factorName = factor(textLabels$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                 "TFIIF", "TFIIE"))


colorPalette = c("#E29930", "#217CA3","#32384D")
#plot
p = ggplot(plotSample, aes(x = time, y = normCount_ends1/Vmax_nls, color = timeGroup))
p  + theme_bw() +
  ylab("HA/Myc") +
  geom_line(aes(x = time, y = y, group = peakName), lwd = 0.8) + 
  geom_point() +
  geom_smooth(aes(x = time, y = ratio, group = peakName),
              method="nls",
              color = "black",
              method.args = list(formula = y ~ 1*x^n / (Jm^n + x^n), 
                                 start=list(Jm=40, n=3)),
              se=FALSE, linetype = "dashed", lwd = 0.8) +
  facet_grid(timeGroup ~ factorName) +
  guides(color = guide_legend(title = paste("residence", "time", sep = "\n")))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size=10),
        strip.background=element_rect(colour="black",
                                      fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = colorPalette) +
  xlab("time [min]") +
  geom_text(data = textLabels, x = 80, y = 0.05, 
            aes(label = myLabels), size = 2.5)






