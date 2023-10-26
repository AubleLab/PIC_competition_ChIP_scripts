rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# get GAL genes
GAL = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(locus=V4, gene = V6)%>% 
  dplyr::filter(startsWith(gene, "GAL")) %>% 
  distinct() %>% 
  dplyr::rename(gene=locus, sysName = gene)


# upload synthesis rates
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  left_join(GAL) %>% 
  mutate(x= "X")  %>% 
  mutate(GAL_class = ifelse(is.na(sysName), "no", "yes"))


synthesisRates$GAL_class = factor(synthesisRates$GAL_class, levels = c("no", "yes"))

p = ggplot(synthesisRates, aes(x = x , y = cropped, fill = GAL_class))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = GAL_class)) + 
  geom_boxplot(outlier.shape=NA, color = "black", alpha = 0.5) +
  theme_classic() +
  xlab(" ") +
  ylab(paste("synthesis rate", "[mRNA/cell/min]", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(size=9),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #facet_wrap(~factorName, scales = "free", nrow = 1) +
  scale_fill_manual(values= c("#E6772E", "#3D4C53")) +
  scale_color_manual(values= c("#E6772E", "#3D4C53")) +
  ylim(0,1.6) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.y = 1.6)

