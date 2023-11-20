#set working directory
setwd("working directory")
library(tidyverse)
library(readxl)

DMPs <- read_csv("Blood_DMPs.csv")
DMPs_hypo <- DMPs %>% filter(Effect < 0) %>% pull(MarkerName)
DMPs_hyper <- DMPs %>% filter(Effect >=  0) %>% pull(MarkerName)
DMPs <- DMPs %>% pull(MarkerName)

VMPs <- read_csv("Blood_VMPs.csv") %>% pull(MarkerName)

DMPs_only <- setdiff(DMPs,VMPs)
DMPs_VMPs <- intersect(DMPs,VMPs)
VMPs_only <- setdiff(VMPs,DMPs)


###############Create background of probes for each clock
all <- read_tsv("Meta-analysis results DMPs blood.txt") %>%
  pull(CpG)
annotation_450K <- read_delim("Annotation_HM450.txt")
annotation_EPIC <- read_delim("Annotation_EPIC.txt")
annotation <- full_join(annotation_450K,
                        annotation_EPIC)

#I have access to Horvath's background, and PhenoAge's background is similar so I'll approximate it as Horvath's
horvath <- read_csv("AdditionalFile22probeAnnotation21kdatMethUsed.csv")$Name
all_horvath <- intersect(all,horvath)

#I have access to Zhang et al., 2019
zhang <- read_delim("blup.coef")$probe
all_zhang <- intersect(all,zhang)

#I have access to the centenarian background CpGs
centenarian <- read_csv("Centenarian clock background CpGs.csv")$CGid
all_centenarian <- intersect(all,centenarian)

#I have access to the universal clock background
universal <- read_delim("GPL28271-57075.txt",
                        comment = "#")$ID
all_universal <- intersect(all,universal)

#I'll approximate the Hannum background as 450K
all_hannum <- intersect(all,annotation_450K$probeID)

#I'll approximate the DunedinPoAm as 450K and EPIC probes
all_poam <- intersect(all,annotation$probeID)

#I'll approximate the DunedinPACE as same as DunedinPoAm but restricted to probes deemed reliable by Sugden et al., 2020 Patterns (Daniel Belsky informed me of it)
reliable_probes <- read_excel("Sugden_MethylationReliability_Data_S1.xlsx",
                              skip = 1) %>%
  filter(Reliability>0.4)%>%
  pull(`Illumina Probe ID`)
all_pace <- intersect(all_poam,reliable_probes)

###############Create chi-square test matrix for each clock  
library(readxl)
sheets <- excel_sheets(path = "Clock sites.xlsx")
sheets <- sheets [-which(sheets%in%c("Universal clock 1",
                                     "Universal clock 3",
                                     "Skin & Blood clock"))]

#### Hypo and hyperDMPs enrichment ####
residuals <- cbind()
pval <- c()
age <- DMPs

for (s in sheets)
{
  clock <- read_excel("Clock sites.xlsx",
                      sheet = s) %>% pull(1)
  
  if (s=="Pan-tissue"|s=="PhenoAge")
  {
    nonclock <- setdiff(all_horvath,clock)
    nonage <- setdiff(all_horvath,age)
      }
  else if (s=="Zhang et al. 2019")
  {
    nonclock <- setdiff(all_zhang,clock)
    nonage <- setdiff(all_zhang,age)
    
  }
  else if (s=="Hannum")
  {
    nonclock <- setdiff(all_hannum,clock)
    nonage <- setdiff(all_hannum,age)
  }
  else if (s=="Centenarian clock")
  {
    nonclock <- setdiff(all_centenarian,clock)
    nonage <- setdiff(all_centenarian,age)
  }
  else if (s=="DunedinPoAm")
  {
    nonclock <- setdiff(all_poam,clock)
    nonage <- setdiff(all_poam,age)
  }
  else if (s=="DunedinPACE")
  {
    nonclock <- setdiff(all_pace,clock)
    nonage <- setdiff(all_pace,age)
  }
  else if (str_detect(s,"Universal"))
  {
    nonclock <- setdiff(all_universal,clock)
    nonage <- setdiff(all_universal,age)
  }
  
 
  #####Chi-square test on unadjusted analyses
  matrix <- cbind(
    #Non-age-related
    c(length(intersect(clock,nonage)),
      length(setdiff(nonage,clock))),
    #HypoDMPs
    c(length(intersect(clock,DMPs_hypo)),
      length(intersect(nonclock,DMPs_hypo))),
    #HyperDMPs
    c(length(intersect(clock,DMPs_hyper)),
      length(intersect(nonclock,DMPs_hyper))))
 
  colnames(matrix) <- c("Non-age-related",
                        "HypoDMPs",
                        "HyperDMPs")
  rownames(matrix) <- c("CpGs in the clock",
                        "CpGs not in the clock")
  chisq <- chisq.test(matrix)
  pval <- c(pval,chisq$p.value)
    residuals <- cbind(residuals,chisq$stdres[1,])
}
colnames(residuals) <- sheets

FDR <- p.adjust(pval)
names(FDR) <- sheets
residuals <- residuals[sort(rownames(residuals)),]
library(corrplot)
library(RColorBrewer)
library(pheatmap)
myBreaks <- c(seq(min(residuals,na.rm=T), 0, length.out=ceiling(100/2) + 1),
              seq(max(residuals,na.rm=T)/100, max(residuals,na.rm=T), length.out=floor(100/2)))

#Add annotation for each clock: clock type, and FDR
clockgen <- data.frame(`Clock type` = c(rep("Chronological age",5),"Biological age",rep("Pace of ageing",2)))
rownames(clockgen) <- sheets
names(clockgen)[names(clockgen)=="Clock.type"] <- "Clock type"

#Add star in each cell if abs(residual) > 2 (indicates a significant depletion or enrichment)
labels <- abs(residuals)
labels[abs(residuals) > 1.96] <- "\u2217"
labels[abs(residuals) < 1.96] <- ""
#Remove * for DunedinPACE
#labels["Constant VMPs","DunedinPACE"] <- ""
#labels["DMPs-VMPs","DunedinPoAm"] <- ""
#labels["Non-age-related","DunedinPoAm"] <- ""
labels["Non-age-related","Universal clock 2"] <- ""
labels["HyperDMPs","Universal clock 2"] <- ""
labels["HypoDMPs","Hannum"] <- ""

tiff("Correlation plot of residuals for clock sites hyper and hypo DMPs.tiff",
     height = 5,
     width = 9,
     res = 300,
     unit = 'in')
pheatmap(residuals,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         display_numbers = labels,
         #gaps_row = c(2,4,6),
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5,6),
         breaks = myBreaks,
         annotation_colors = list(
           `Clock type` = c(
             "Chronological age" = "#8A66D9",
             "Biological age" = "#BAA0F2",
             "Pace of ageing" = "#CEBDF2"
           )),
         annotation_col = clockgen,
         show_rownames = T)
dev.off()



#### homoDMP, DMP-VMP, contVMP enrichment ####
residuals <- cbind()
pval <- c()
age <- union(DMPs,VMPs)

for (s in sheets)
{
  clock <- read_excel("Clock sites.xlsx",
                      sheet = s) %>% pull(1)
  
  if (s=="Pan-tissue"|s=="PhenoAge")
  {
    nonclock <- setdiff(all_horvath,clock)
    nonage <- setdiff(all_horvath,age)
    #nonage_CTC <- setdiff(all_horvath,age_CTC)
  }
  else if (s=="Zhang et al. 2019")
  {
    nonclock <- setdiff(all_zhang,clock)
    nonage <- setdiff(all_zhang,age)
    #nonage_CTC <- setdiff(all_zhang,age_CTC)
  }
  else if (s=="Hannum")
  {
    nonclock <- setdiff(all_hannum,clock)
    nonage <- setdiff(all_hannum,age)
    #nonage_CTC <- setdiff(all_hannum,age_CTC)
  }
  else if (s=="Centenarian clock")
  {
    nonclock <- setdiff(all_centenarian,clock)
    nonage <- setdiff(all_centenarian,age)
    #nonage_CTC <- setdiff(all_centenarian,age_CTC)
  }
  else if (s=="DunedinPoAm")
  {
    nonclock <- setdiff(all_poam,clock)
    nonage <- setdiff(all_poam,age)
    #nonage_CTC <- setdiff(all_poam,age_CTC)
  }
  else if (s=="DunedinPACE")
  {
    nonclock <- setdiff(all_pace,clock)
    nonage <- setdiff(all_pace,age)
    #nonage_CTC <- setdiff(all_pace,age_CTC)
  }
  else if (str_detect(s,"Universal"))
  {
    nonclock <- setdiff(all_universal,clock)
    nonage <- setdiff(all_universal,age)
    #nonage_CTC <- setdiff(all_universal,age_CTC)
  }
  
  #####Chi-square test on unadjusted analyses
  matrix <- cbind(
    #Non-age-related
    c(length(intersect(clock,nonage)),
      length(setdiff(nonage,clock))),
    #DMPs only
    c(length(intersect(clock,DMPs_only)),
      length(intersect(nonclock,DMPs_only))),
    #VMPs only
    c(length(intersect(clock,VMPs_only)),
      length(intersect(nonclock,VMPs_only))),
    #DMPs/VMPs
    c(length(intersect(clock,DMPs_VMPs)),
      length(intersect(nonclock,DMPs_VMPs))))
  
  colnames(matrix) <- c("Non-age-related",
                        "Homoscedastic DMPs",
                        "Constant VMPs",
                        "DMPs-VMPs")
  rownames(matrix) <- c("CpGs in the clock",
                        "CpGs not in the clock")
  chisq <- chisq.test(matrix)
  pval <- c(pval,chisq$p.value)
  
  residuals <- cbind(residuals,chisq$stdres[1,])
}
colnames(residuals) <- sheets

FDR <- p.adjust(pval)
names(FDR) <- sheets
# names(FDR) <- rep(sheets,
#                   each = 2)

#Remove non-age-related?
#residuals <- residuals[-grep("Non",rownames(residuals)),]
residuals <- residuals[sort(rownames(residuals)),]
library(corrplot)
library(RColorBrewer)
library(pheatmap)
myBreaks <- c(seq(min(residuals,na.rm=T), 0, length.out=ceiling(100/2) + 1),
              seq(max(residuals,na.rm=T)/100, max(residuals,na.rm=T), length.out=floor(100/2)))

#Add annotation for each clock: clock type, and FDR
clockgen <- data.frame(`Clock type` = c(rep("Chronological age",5),"Biological age",rep("Pace of ageing",2)))
rownames(clockgen) <- sheets
names(clockgen)[names(clockgen)=="Clock.type"] <- "Clock type"

#Add star in each cell if abs(residual) > 2 (indicates a significant depletion or enrichment)
labels <- abs(residuals)
labels[abs(residuals) > 1.96] <- "\u2217"
labels[abs(residuals) < 1.96] <- ""
#Remove * for DunedinPACE
labels["Constant VMPs","DunedinPACE"] <- ""
labels["DMPs-VMPs","DunedinPoAm"] <- ""
labels["Non-age-related","DunedinPoAm"] <- ""


tiff("Correlation plot of residuals for clock sites hypo DMP, DMP-VMP, const VMPs.tiff",
     height = 5,
     width = 9,
     res = 300,
     unit = 'in')
pheatmap(residuals,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         display_numbers = labels,
         #gaps_row = c(2,4,6),
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5,6),
         breaks = myBreaks,
         annotation_colors = list(
           `Clock type` = c(
             "Chronological age" = "#8A66D9",
             "Biological age" = "#BAA0F2",
             "Pace of ageing" = "#CEBDF2"
           )),
         annotation_col = clockgen,
         show_rownames = T)
dev.off()



