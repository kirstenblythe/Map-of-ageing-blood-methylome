#set working directory
setwd("working directory")
library(tidyverse)


#Upload list of DMPs and VMPs from meta-analysis
DMPs <- read_csv("Blood_DMPs.csv")
VMPs <- read_csv("Blood_VMPs.csv") 

#### hyperDMP and hypoDMP chromatin state enrichment ####
#classify hypo and hyper DMPs

all <- read_tsv("Meta-analysis results DMPs blood.txt") %>%
  pull(CpG)

DMPs_hypo <- DMPs %>% 
  filter(Effect <0) %>% 
  pull(MarkerName)
DMPs_hyper <- DMPs %>% 
  filter(Effect >=0 ) %>% 
  pull(MarkerName)

DMPs <- DMPs %>% pull(MarkerName)

#remaining non-age-related CpGs
nonage <- setdiff(all,c(DMPs_hypo,DMPs_hyper))

#Create tables to merge with annotation
tomerge <- tibble(probeID = c(DMPs_hypo,DMPs_hyper, nonage),
                  class = c(rep("hypoDMP",length(DMPs_hypo)),
                            rep("hyperDMP",length(DMPs_hyper)),
                            rep("Nonage",length(nonage))))

#Load annotation
annotation <- read_delim("AnnotationBlood.txt") %>%
  select(probeID,E062) %>%
  unique() %>%
  filter(str_detect(E062, #chromatin states in PBMCs
                    ",",
                    negate = T))

annotation$E062 <- gsub('[[:digit:]]+',"",annotation$E062)
annotation$E062 <- gsub('_',"",annotation$E062)
annotation$E062[annotation$E062 == "TssA"] <- "Active TSS"
annotation$E062[annotation$E062 == "TssAFlnk"] <- "Flanking active TSS"
annotation$E062[annotation$E062 == "TxFlnk"] <- "Transcr. at gene 5' and 3'"
annotation$E062[annotation$E062 == "Tx"] <- "Strong transcription"
annotation$E062[annotation$E062 == "TxWk"] <- "Weak transcription"
annotation$E062[annotation$E062 == "EnhG"] <- "Genic enhancers"
annotation$E062[annotation$E062 == "Enh"] <- "Enhancers"
annotation$E062[annotation$E062 == "ZNF/Rpts"] <- "ZNF genes & repeats"
annotation$E062[annotation$E062 == "Het"] <- "Heterochromatin"
annotation$E062[annotation$E062 == "TssBiv"] <- "Bivalent/poised TSS"
annotation$E062[annotation$E062 == "BivFlnk"] <- "Flanking bivalent TSS/Enh"
annotation$E062[annotation$E062 == "EnhBiv"] <- "Bivalent enhancer"
annotation$E062[annotation$E062 == "ReprPC"] <- "Repressed polycomb"
annotation$E062[annotation$E062 == "ReprPCWk"] <- "Weak repressed polycomb"
annotation$E062[annotation$E062 == "Quies"] <- "Quiescent/low"

chrom_state <- inner_join(tomerge,
                          annotation) %>%
  group_by(E062,
           class) %>%
  tally() %>%
  mutate(perc = case_when(class == "hypoDMP" ~ round((n/length(DMPs_hypo))*100,1),
                          class == "hyperDMP" ~ round((n/length(DMPs_hyper))*100,1),
                          class == "Nonage" ~ round((n/length(nonage))*100,1)))

tiff("Distribution of hyper and hypo DMPs in chromatin states.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300)
ggplot(data = chrom_state) +
  geom_col(mapping = aes(x = E062,
                         y = perc,
                         fill = class), 
           position = "dodge") +
  xlab("") + 
  scale_fill_manual(values = c("#63C3A1","#A3D8C3","#D9D9D9"),
                    name = "CpG classification",
                    labels = c("hyper", "hypo", "non-DMP")) +
  #coord_flip() +
  #scale_x_discrete(limits=rev) +
  theme_minimal()+
  theme(# Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 7.5, angle = 90, hjust = 0.5),
    axis.text.y = element_text(color = "gray12", size = 7.5),
    panel.grid = element_blank()
  ) +
  labs(title = "Distribution of DMPs in chromatin states")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.y = element_text(size = 10, vjust = 1.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))
dev.off()



####Residuals
forx2_chrstate <- chrom_state %>% 
  dplyr::select(E062, class, n) %>% 
  drop_na(E062) %>% 
  pivot_wider(id_cols = E062,
              names_from = class,
              values_from = n)
class(forx2_chrstate) 
forx2_chrstate <- as.data.frame(forx2_chrstate)
rownames(forx2_chrstate) <- forx2_chrstate$E062
forx2_chrstate <- forx2_chrstate %>% dplyr::select(-E062)

forx2_chrstate <- t(forx2_chrstate)
chisq_chrstate <- chisq.test(forx2_chrstate)


library(corrplot)
tiff("Correlation plot of x2 residuals for chromatin states of the three classes of age-related CpGs.tiff",
     height = 5,
     width = 8,
     res = 300,
     unit = 'in')
corrplot(chisq_chrstate$residuals, 
         method = "square",
         is.cor = FALSE, 
         cl.ratio=1,
         tl.col="black",
         cl.pos = "n",
         col = colorRampPalette(c("#5EA7E6","white","#5758F2"))(100),
         tl.cex = .75)
dev.off()


#### homoscedastic DMP, DMP-VMP and constant VMP chromatin state enrichment ####
#upload entire list of CpGs meta-analysed in DMP meta

all <- read.table("Meta-analysis results DMPs blood.txt") %>%
  pull(CpG)

DMPs <- DMPs %>% pull(MarkerName)
VMPs <- VMPs %>% pull(MarkerName)

#separate into three groups: homoscedastic DMPs, DMPs-VMPs, constant VMPs

#homoscedastic DMPs
DMPs_only <- setdiff(DMPs,VMPs)

#DMPs-VMPs
DMPs_VMPs <- intersect(DMPs,VMPs)

#constant VMPs
VMPs_only <- setdiff(VMPs,DMPs)

#remaining non-age-related CpGs
nonage <- setdiff(all,c(DMPs_only,DMPs_VMPs,VMPs_only))

#Create tables to merge with annotation
tomerge <- tibble(probeID = c(DMPs_only,DMPs_VMPs,VMPs_only,nonage),
                  class = c(rep("HomoDMP",length(DMPs_only)),
                            rep("DMP-VMP",length(DMPs_VMPs)),
                            rep("ConstVMP",length(VMPs_only)),
                            rep("Nonage",length(nonage))))

chrom_state <- inner_join(tomerge,
                          annotation) %>%
  group_by(E062,
           class) %>%
  tally() %>%
  mutate(perc = case_when(class == "HomoDMP" ~ round((n/length(DMPs_only))*100,1),
                          class == "DMP-VMP" ~ round((n/length(DMPs_VMPs))*100,1),
                          class == "ConstVMP" ~ round((n/length(VMPs_only))*100,1),
                          class == "Nonage" ~ round((n/length(nonage))*100,1)))

tiff("Distribution of the three classes of age-related CpGs in PBMC chromatin states.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300)
ggplot(data = chrom_state) +
  geom_col(mapping = aes(x = E062,
                         y = perc,
                         fill = class), 
           position = "dodge")+
  ylab("% of CpGs") +
  xlab("") + 
  scale_fill_manual(values = c("#2d2e74","#8250a0","#81bfe9","#D9D9D9"),
                    name = "CpG classification",
                    labels = c("ConstVMP", "DMP-VMP", "HomoDMP","Non-age")) +
  theme_minimal()+
  theme(# Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 7.5, angle = 90, hjust = 0.5),
    axis.text.y = element_text(color = "gray12", size = 7.5),
    panel.grid = element_blank()
  ) +
  labs(title = "Distribution of age-related CpGs in PBMC chromatin states")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.y = element_text(size = 10, vjust = 1.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))
dev.off()


####Residuals
forx2_chrstate <- chrom_state %>% 
  dplyr::select(E062, class, n) %>% 
  drop_na(E062) %>% 
  pivot_wider(id_cols = E062,
              names_from = class,
              values_from = n)
class(forx2_chrstate) 
forx2_chrstate <- as.data.frame(forx2_chrstate)
rownames(forx2_chrstate) <- forx2_chrstate$E062
forx2_chrstate <- forx2_chrstate %>% dplyr::select(-E062)

forx2_chrstate <- t(forx2_chrstate)
chisq_chrstate <- chisq.test(forx2_chrstate)


library(corrplot)
tiff("Correlation plot of x2 residuals for chromatin states of the three classes of age-related CpGs.tiff",
     height = 5,
     width = 8,
     res = 300,
     unit = 'in')
corrplot(chisq_chrstate$residuals, 
         method = "square",
         is.cor = FALSE, 
         cl.ratio=1,
         tl.col="black",
         cl.pos = "n",
         col = colorRampPalette(c("#457ABF","white","#914BBF"))(100),
         tl.cex = .75)
dev.off()

