rm(list = ls())

library(tidyverse)
library(viridis)

# upload table with residence times replace <1 min with a random number between 0-1
# make class labels
resTimes = read.csv("data/residence_times_all.csv") %>% 
  dplyr::select(-peakName) %>% 
  pivot_longer(cols = TBP:TFIIF, names_to = "factorName", values_to = "resTime", values_drop_na = TRUE) %>% 
  distinct()

countTimeGroups = resTimes %>% 
  mutate(randNum = runif(nrow(resTimes), min=0, max=1)) %>% 
  mutate(resTimesNum = ifelse(resTime == "<1", randNum, as.numeric(resTime)))%>% 
  mutate(timeGroup = ifelse(resTimesNum < 1, "<1 min",
                                   ifelse(resTimesNum < 10, "1-10 min", ">10 min"))) %>% 
  dplyr::select(factorName, timeGroup) %>% 
  group_by(factorName, timeGroup) %>% 
  tally()


countTimeGroups$factorName = factor(countTimeGroups$factorName, levels = c("TBP", "TFIIA", "TFIIB", 
                                                                       "TFIIF", "TFIIE"))
countTimeGroups$timeGroup = factor(countTimeGroups$timeGroup, 
                                   levels = c("<1 min", "1-10 min", ">10 min"))
colorPalette = c("#E29930", "#217CA3","#32384D")
p = ggplot(countTimeGroups, aes(x = factorName, y = n, fill = timeGroup))
p + geom_bar(stat= "identity") +
  scale_fill_manual(values = colorPalette) +
  theme_classic() +
  xlab("")+
  guides(fill = guide_legend(title = paste("residence", "time", sep = "\n")))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(face = "bold", hjust = 1,
                                   size = 9, angle = 45),
        text = element_text(size=10)) 

