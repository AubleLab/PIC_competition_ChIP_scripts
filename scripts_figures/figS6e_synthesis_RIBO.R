rm(list = ls())

library(tidyverse)
library(ggpubr)

set.seed(42)

# upload matching systematic names from www.yeastgenome.org:
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/dbxref.tab
geneConvertTable = read.delim("data/dbxref.tab", header = F) %>% 
  dplyr::select(V4, V6) %>% 
  dplyr::rename(gene = V4) %>% 
  dplyr::rename(geneSys = V6)


# upload table with synthesis rates and identify ribosomal subunits
synthesisRates = read.csv("data/synthesisRates.csv") %>% 
  mutate(cropped = ifelse(synthesis_perCell_perMin > 1.5, 1.5, synthesis_perCell_perMin))%>% 
  left_join(geneConvertTable) %>% 
  mutate(ribo = ifelse(str_detect(gene, "RPL") | str_detect(geneSys, "RPL"), "yes", "no")) %>% 
  replace_na(list(ribo = "no")) %>% 
  distinct() %>% 
  mutate(x= "X")


synthesisRates$ribo = factor(synthesisRates$ribo, levels = c("no", "yes"))
# plot
p = ggplot(synthesisRates, aes(x = x , y = cropped, fill = ribo))
p + geom_point(position=position_jitterdodge(jitter.width = 0.3), 
               aes(color = ribo), alpha = 0.5) + 
  geom_boxplot(outlier.shape=NA, color = "black") +
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
  scale_fill_manual(values= c("#4DB3B3", "#E64A45")) +
  scale_color_manual(values= c("#4DB3B3", "#E64A45")) +
  ylim(0,1.6) + 
  stat_compare_means(method = "t.test", label = "p.signif", label.y = 1.6)
#ggsave("figures/panels/figS6/synthesis_vs_RIBO.pdf", width = 2.8, height = 6, units = "cm")


