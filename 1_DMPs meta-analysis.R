#### Perform EWAS in each blood dataset ####
#using GSE197676 as example dataset
#dataset has been preprocessed according to ChAMP pipeline

#set directory
path = 'pvol/Preprocessing/Blood/GEO/GSE197676'
setwd(path)

library(tidyverse)
library(limma)

#Identify DMPs in each dataset

create_summary <- function(toptable = NULL,
                           dataset_label = NULL,
                           directory = getwd())
{
  CPG <- rownames(toptable)
  ALLELE1 <- rep(1,nrow(toptable))
  ALLELE2 <- rep(2,nrow(toptable))
  TESTSTAT <- toptable$t
  PVALUE <- toptable$P.Value
  EFFECTSIZE <- toptable$logFC
  SE <- toptable$SE
  results = data.frame(CPG,
                       ALLELE1,l
                       ALLELE2,
                       TESTSTAT,
                       PVALUE,
                       EFFECTSIZE,
                       SE)
  write.table(results,
              file=paste0(directory,"/",dataset_label,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

#load the data

B <- B <- data.table::fread("GSE197676 beta after normalisation and batch correction.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B %>% select(-V1)
M <- logit2(B)
pheno <- read.delim("GSE197676 Phenotypes.txt")
glimpse(pheno)

#run linear model including the dataset specific covariates (Supplementary Table 1)
design=model.matrix(~age +
                      sex,
                    pheno)

#first run the linear model on the M values
fit1_M <- lmFit(M,
                design)
fit2_M <- eBayes(fit1_M)

#repeat the linear model on the beta values to extract the effect size in beta
fit1_B <- lmFit(B,
                design)
fit2_B <- eBayes(fit1_B)

#extract the results for age
coef = "age"
results <- topTable(fit2_M,
                    coef=coef,
                    number = Inf,
                    p.value = 1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"] #effect size in beta values
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#save p-value histogram
#Save p-value histogram
tiff('GSE197676_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue",bins = 30)
dev.off()

#check the number of sig DMPs in the dataset
fdr = 0.005
results_age=topTable(fit2_M,
                     coef = "age",
                     number = nrow(M),
                     adjust.method = "BH",
                     p.value = fdr)
#save the results as a table to run the meta-analysis
directory = "/pvol/EWAS/Blood/DMPs"
create_summary(toptable = results,
               dataset_label = "GSE197676",
               directory = directory)

#save the residuals to run the Breusch-Pagan test for VMPs
resid <- residuals(fit2_M, M)
write.table(signif(resid,digits = 4),
            file="GSE197676_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

#check to see if limma worked


cpg <-rownames(results)[1]
pheno_with_meth <- cbind(pheno, meth= as.numeric(B[cpg,]))
tiff('GSE197676_DMP_check.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(pheno_with_meth, aes(x=age, y=meth)) +
  geom_jitter(width = 0.06)+
  labs(y=paste("Methylation at",cpg))
dev.off()

#repeat for remaining datasets 

#### formatting files for METAL ####

#### This is just to explain how the files were formatted before the meta-analysis was run. We cannot provide the raw data for each dataset so skip this step 
library(bacon)
library(tidyverse)

#set working directory
setwd("/pvol/EWAS/Blood/DMPs")

#List all .tbl files 
files <- list.files()[grep(".tbl",list.files())]
files 

#Run each file in bacon using a loop
for (f in files)
{
    file <- read_tsv(f) #Read the file
    f <- sub("\\.tbl", "", f) #Obtain the name of the dataset without the ".tbl" at the end
    print(f) #show which dataset we are currently analysing
    bc <- bacon(teststatistics = NULL, #run bacon on effect sizes and standard errors
                effectsizes = file$EFFECTSIZE,
                standarderrors = file$SE)
    
    tiff(paste0('QQ-plot_',f,'.tiff'), #save the graph of raw and adjusted effect sizes and standard errors
         width =4,
         height = 2.5,
         units = 'in',
         res=600)
    print(plot(bc, type="qq")) #q-q plot shows the deviation of p-value distribution from null hypothesis. The observed P values for each CpG are sorted from largest to smallest and plotted against expected values from a theoretical ??2-distribution.
    dev.off()
    
    print(c(inflation(bc), #show the inflation factor for this dataset
            bias(bc))) #show the bias factor for this dataset
    
    file$EFFECTSIZE_CORR <- es(bc)[,1] #add a column to the file corresponding to the corrected effect size
    file$SE_CORR <- se(bc)[,1] #add a column to the file corresponding to the corrected standard error
    file$PVALUE_CORR <- pval(bc)[,1] #add a column to the file corresponding to the corrected p-value
    write.table(file, #save the file
                file=paste0(getwd(),"/",f,"_corrected.tbl"),
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                sep="\t")
}

#multiply the effect size and coefficient for each file because numbers are small
files <- list.files()[grep("_corrected.tbl",list.files())]
files

library(tidyverse)
for (f in files)
{
  file <- read_tsv(f) #Read the file
  f <- sub("\\.tbl", "", f) #Obtain the name of the dataset without the ".tbl" at the end
  print(f) #show which dataset we are currently analysing
  
  file$EFFECTSIZE_CORR <- file$EFFECTSIZE_CORR*100 #add a column to the file corresponding to the corrected effect size
  file$SE_CORR <- file$SE_CORR*100 #add a column to the file corresponding to the corrected standard error
  write.table(file, #save the file
              file=paste0(getwd(),"/",f,"2.tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

#Write code for METAL
#List tbl files
files <- list.files(pattern="2.tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME STDERR",
                  "MARKER CPG",
                  "EFFECT EFFECTSIZE_CORR",
                  "SEPARATOR TAB",
                  "STDERR SE_CORR",
                  "PVALUE PVALUE_CORR",
                  sep = "\r"),
            file = file,
            quote=F,
            row.names=F,
            col.names=F)
#Add datasets to analyse
add_datasets_METAl <- function(d = NULL,
                               f = NULL)
{
  write(paste("PROCESS",
              d,
              sep = "\t"),
        file=f,
        append=TRUE)
}
sapply(files,
       add_datasets_METAl,
       f = file)

#Add last line of code
write("ANALYZE HETEROGENEITY",
      file=file,
      append=TRUE)



#### BEGIN ANALYSIS FROM HERE ####

#### METAL results ####

library(tidyverse)
setwd("working directory")

list.files()[grep("METAANALYSIS",list.files())]

Meta_DMP <- read.table("METAANALYSIS1.TBL",
                   header = TRUE)

#filter on CpGs that are present in at least 3 datasets 
Meta_DMP <- Meta_DMP %>% mutate("Number of studies" = HetDf + 1)
Meta_DMP <- Meta_DMP %>% filter(`Number of studies` >= 3)

#pvalue adjustment
Meta_DMP <- Meta_DMP %>% mutate(pval_adj = p.adjust(Meta_DMP$P.value, method = "fdr", n = length(Meta_DMP$P.value)))

#Filter FDR < 0.005
DMPs <- filter(Meta_DMP, pval_adj < 0.005)

write.csv(DMPs, "Blood_DMPs.csv")

#hypermethylated
sum(sign(DMPs$Effect >= 0))
#113047 hyperDMPs

#hypomethylated
sum(sign(DMPs$Effect < 0))
#220253

#visualise results using volcano plot 
Meta_DMP <- Meta_DMP %>% 
  mutate(Sig = ifelse(pval_adj < 0.005, "FDR < 0.005","Not sig"))
Meta_DMP <- Meta_DMP %>% 
  mutate(coef = ifelse(Effect >= 0, "Positive", "Negative"))
Meta_DMP = Meta_DMP %>% mutate(Classification=ifelse(Sig == "Not sig","nonDMP",
                                                     ifelse(coef=="Negative","hypoDMP","hyperDMP")))

#assign the pvalues that are 0 the lowest possible value in R for volcano plot 

Meta_DMP <- Meta_DMP %>% mutate(P.value_volcano = ifelse(P.value == "0", "5e-324", P.value))
Meta_DMP$P.value_volcano <- as.numeric(Meta_DMP$P.value_volcano)

tiff('Blood DMPs volcano plot.tiff',
     width = 8,
     height = 9,
     units = 'in',
     res = 600)
ggplot(Meta_DMP, aes(Effect, -log10(P.value_volcano))) +
  geom_point(aes(col = Classification), size = 0.05) +
  scale_color_manual(values=c("#63C3A1", "#A3D8C3","grey"),
                     name = "") +
  xlab("% DNA methylation change per year") +
  ylab("-log10(P.value)") +
  theme_minimal() +
  theme(legend.key.height = unit(2,"cm"),
        legend.key.width = unit(3,"cm"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
dev.off()

#### Forest plots ####


#Create empty list
L.names <- Meta_DMP$MarkerName
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(CpG = Meta_DMP$MarkerName,
              Dataset = rep("Meta-analysis",nrow(Meta_DMP)),
              ES = Meta_DMP$`Effect`,
              SE = Meta_DMP$`StdErr`,
              PVAL = Meta_DMP$`P.value`,
              FDR = Meta_DMP$pval_adj,
              Type = rep("Meta-analysis",nrow(Meta_DMP)))
tib <- tib %>%
  mutate_if(is.numeric,
            signif,
            digits = 4)
#Initiate the list from meta-analysis
L <- split(tib, seq(nrow(tib))) #splits each row of the tibble into a list component
names(L) <- Meta_DMP$MarkerName

#load the names of the datasets

library(readxl)
list.files()
datasets <- read_excel("/pvol/EWAS/Blood/Meta-analysis/Dataset names for meta-analysis.xlsx")

#Add each study to each component of the list

setwd("/pvol/EWAS/Blood/DMPs")

for (s in datasets$Dataset)
{
  print(s)
  file <- read_tsv(paste0(s,"_corrected2.tbl"))
  
  #Add FDR to the table
  file <- file %>%
    mutate(FDR = p.adjust(PVALUE))
  
  #Add Effect 
  #Create tibble
  summary <- tibble(CpG = file$CPG,
                    Dataset = rep(s,nrow(file)),
                    ES = file$EFFECTSIZE_CORR,
                    SE = file$SE_CORR,
                    PVAL = file$PVALUE_CORR,
                    FDR = file$FDR,
                    Type = rep("Individual study",nrow(file)))
  summary <- summary %>%
    mutate_if(is.numeric,
              signif,
              digits = 2)
  
  summary <- summary %>%
    filter(CpG %in% names(L))
  subL <- split(summary, seq(nrow(summary)))
  names(subL) = summary$CpG
  
  #Merge the pieces of the two lists that are in common
  L2 <- L[names(subL)]
  L3 <- L[setdiff(names(L),names(subL))]
  listinter <- Map(bind_rows,
                   L2,
                   subL)
  L <- c(listinter,
         L3)
}

saveRDS(L,
       file = "List for meta-analysis forest plot blood DMPs.rds")

 
#visualise most significant DMP as a forest plot
L <- readRDS("List for meta-analysis forest plot blood DMPs.rds")

cg <- L[[1]]
cg$CpG

Datasets <- datasets %>% select(Dataset,N,Array)

cg <- left_join(cg, Datasets, by = "Dataset")

cg$N <- as.numeric(cg$N)
N_meta <- sum(cg$N, na.rm = TRUE)
cg$N[is.na(cg$N)]<-sum(cg$N, na.rm = T)

cg_Meta <- filter(cg, Dataset == "Meta-analysis")
cg <- cg %>% filter(Dataset != "Meta-analysis")

tiff("cg16867657_forest.tiff",
     height = 11,
     width = 7,
     unit = 'in',
     res = 500)
forest(cg$ES,
       sei = cg$SE,
       slab = cg$Dataset,
       xlab = "% DNAm change per year",
       psize = 1,
       header = "Dataset",
       ilab = cg$N,
       ilab.xpos = -0.1,
       order = -(cg$N),
       cex = 0.6,
       xlim = c(-0.35,1.05))
addpoly(cg_Meta$ES,
        sei = cg_Meta$SE,
        cex = 0.6,
        mlab = "Meta-analysis",
        col = "Purple",
        border = "Purple")
text(-0.1, 53,
     "N",
     cex = 0.8)
dev.off()

#### Methylation plots ####

#Visualise most significant cpg across three independent datasets 

#JHS

#the beta file and phenotype file have been provided 

JHS_B <- data.table::fread("/pvol/Preprocessing/Blood/JHS/JHS beta.txt")
JHS_B <- as.data.frame(JHS_B)
rownames(JHS_B) <- JHS_B$V1
JHS_B <- JHS_B %>% select(-V1)
JHS_pheno <- read.delim("/pvol/Preprocessing/Blood/JHS/JHS Phenotypes.txt", sep = "")

cpg <-"cg16867657"
pheno_with_meth <- cbind(JHS_pheno, meth= as.numeric(JHS_B[cpg,]))
tiff('JHS_cg16867657_grey.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 400)
ggplot(pheno_with_meth, aes(x=age, y=meth*100)) +
  geom_jitter(width = 0.2, colour = "#63C3A1")+
  labs(y=paste("MF %")) +
  geom_smooth(method = "lm", colour = "grey") +
  xlim(19,106) + 
  ggtitle("JHS") +
  theme_minimal() 
dev.off()

#GSE40279

GSE40279_B <- data.table::fread("/pvol/Preprocessing/Blood/GEO/GSE40279/GSE40279 beta after filtering and imputation.txt")
GSE40279_B <- as.data.frame(GSE40279_B)
rownames(GSE40279_B) <- GSE40279_B$V1
GSE40279_B <- GSE40279_B %>% select(-V1)
GSE40279_pheno <- read.csv("/pvol/Preprocessing/Blood/GEO/GSE40279/GSE40279 pheno with cell types.csv")
cpg <-"cg16867657"
pheno_with_meth <- cbind(GSE40279_pheno, meth= as.numeric(GSE40279_B[cpg,]))
tiff('GSE40279_cg16867657.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 300)
ggplot(pheno_with_meth, aes(x=age, y=meth*100)) +
  geom_jitter(width = 0.2, colour = "#63C3A1")+
  labs(y=paste("MF %")) +
  geom_smooth(method = "lm", colour = "grey") +
  xlim(19,106) + 
  ggtitle("GSE40279") +
  theme_minimal() +
  theme(title = element_text(hjust = 0.5))
dev.off()


#GSE152026

GSE152026_B <- data.table::fread("/pvol/Preprocessing/Blood/GEO/GSE152026/GSE152026 beta normalised and batch corrected.txt")
GSE152026_B <- as.data.frame(GSE152026_B)
rownames(GSE152026_B) <- GSE152026_B$V1
GSE152026_B <- GSE152026_B %>% select(-V1)
GSE152026_pheno <- read.delim("/pvol/Preprocessing/Blood/GEO/GSE152026/GSE152026 Phenotypes.txt", sep = "")
GSE152026_B <- GSE152026_B[,GSE152026_pheno$Sample_Name]
cpg <-"cg16867657"
pheno_with_meth <- cbind(GSE152026_pheno, meth= as.numeric(GSE152026_B[cpg,]))
tiff('GSE152026_cg16867657.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 600)
ggplot(pheno_with_meth, aes(x=age, y=meth*100)) +
  geom_jitter(width = 0.2, colour = "#63C3A1")+
  labs(y=paste("MF %")) +
  geom_smooth(method = "lm", colour = "grey") +
  xlim(19,106) + 
  ggtitle("GSE152026") +
  theme_minimal() 
dev.off()


#### Annotation ####

#load annotation

annotation <- read.delim("AnnotationBlood.txt")
#Change name to probeID
Meta_DMP <- dplyr::rename(Meta_DMP,
                          probeID = MarkerName)

meta_anno <- left_join(Meta_DMP, annotation, by= "probeID")

meta_anno <- meta_anno %>% 
  dplyr::select(probeID,
                CpG_chrm,
                CpG_beg,
                CGIposition,
                E062,
                genesUniq_with_enh,
                Effect,
                StdErr,
                `P.value`,
                pval_adj,
                HetDf,
                HetISq,
                HetPVal,
                Sig,
                coef,
                Classification,
                P.value_volcano)

meta_anno <- meta_anno %>% 
  mutate(genesUniq_with_enh = replace_na(genesUniq_with_enh,""))

#Change HetDf by number of included studies
meta_anno$HetDf <- meta_anno$HetDf +1

colnames(meta_anno) = c("CpG",
                         "Chromosome",
                         "Position (hg38)",
                         "CpG island position",
                         "Chromatin state in primary mononuclear cells from peripheral blood",
                         "Annotated gene(s)",
                         "Effect size",
                         "SE",
                         "P-value",
                         "FDR",
                         "Number of studies",
                         "Heterogeneity index (I2)",
                         "Heterogeneity p-value",
                         "Significance",
                         "Direction",
                         "Classification",
                         "P-value volcano")

write.table(meta_anno,
            file="1_Meta-analysis results DMPs blood.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#### Classification of entropic vs antientropic DMPs ####
library(tidyverse)

BIOS_B <- data.table::fread("/pvol/Preprocessing/Blood/BIOS/BIOS Normalized beta values after batch correction.txt")
BIOS_B <- as.data.frame(BIOS_B)
rownames(BIOS_B) <- BIOS_B$V1
BIOS_B <- BIOS_B %>% select(-V1)
BIOS_pheno <- read.delim("/pvol/Preprocessing/Blood/BIOS/BIOS Phenotypes.txt", sep = "")
meta_total <- read.delim("/pvol/EWAS/Blood/Meta-analysis/DMPs/Meta-analysis results DMPs blood.txt",
                         sep = "")

#filter the young samples to calculate methylation levels at baseline
young_samples <- BIOS_pheno %>% filter(age < 30)
#475 young samples with age < 30

#Calculate the avergae methylation fraction (MF%) at 'baseline' 
BIOS_B_young_avg <- as.data.frame(rowMeans(BIOS_B[,young_samples$Sample_Name]))
names(BIOS_B_young_avg)[names(BIOS_B_young_avg)=="rowMeans(BIOS_B[, young_samples$Sample_Name])"] <- "avg MF%"
BIOS_B_young_avg$CpG <- rownames(BIOS_B_young_avg)
BIOS_B_young_avg <- as_tibble(BIOS_B_young_avg) %>% left_join(meta_total, by = "CpG")

#classify meth level
BIOS_meta_allcpgs <- BIOS_B_young_avg %>% mutate("meth level" = case_when(`avg MF%` >= 0.75 ~ "High",
                                                                          `avg MF%` <= 0.25 ~ "Low",
                                                                          `avg MF%` < 0.75 | `avg MF%` > 0.25 ~ "Intermediate"))
#percent of meth high, low or intermed methylated 
BIOS_meta_allcpgs %>% group_by(`meth level`) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))

#46% highly methylated
#40% lowly methylated
#13% intermediate methylation 

#filter on DMPs
BIOS_meta_DMPs <- BIOS_meta_allcpgs %>% filter(FDR < 0.005)

#percent of meth high, low or intermed methylated 
BIOS_meta_DMPs %>% group_by(`meth level`) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))

#47% of DMPs highly methylated
#32% lowly methylated 
#21% intermediate methylated

#what proportion of highly methylated CpGs lose methylation with age?
BIOS_meta_DMPs %>% filter(`meth level` == "High") %>% group_by(Direction) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
#roughly 73% of highly methylated CpGs in young lose methylation with age and the remaining 23% trend towards 100% methylation

#what proportion of intermediate methylated CpGs gain methylation with age?
BIOS_meta_DMPs %>% filter(`meth level` == "Intermediate") %>% group_by(Direction) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
#40% gain methylation 
#60% lose methylation
#what proportion of intermediate DMPs are >50% 
BIOS_meta_DMPs %>% filter(`meth level` == "Intermediate") %>% mutate(above_50 = case_when(`avg MF%` >0.5 ~ "51-74%",
                                                                                          `avg MF%` <0.5 ~ "26-49%")) %>% group_by(above_50) %>% summarise(n = n()) %>% mutate(freq = n / sum(n)) 
#45% 26-49% MF
#55% 51-74% MF

#what proportion of lowly methylated CpGs gain methylation with age?
BIOS_meta_DMPs %>% filter(`meth level` == "Low") %>% group_by(Direction) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
#82% gain methylation

#add the avg MF% of the old samples to determine if the avg MF of the intermed samples increases or decreases with age 
old_samples <- BIOS_pheno %>% filter(age >60)
#104 young samples with age > 60

#calculate avg MF% for old
BIOS_B_old_avg <- as.data.frame(rowMeans(BIOS_B[,old_samples$Sample_Name]))
names(BIOS_B_old_avg)[names(BIOS_B_old_avg)=="rowMeans(BIOS_B[, old_samples$Sample_Name])"] <- "avg MF% old"
BIOS_B_old_avg$CpG <- rownames(BIOS_B_old_avg)
BIOS_meta_DMPs <- left_join(BIOS_meta_DMPs, BIOS_B_old_avg, by = "CpG")

#classify the DMPs according to whether the DMPs gained or lost meth during ageing 
BIOS_meta_DMPs <- BIOS_meta_DMPs %>% mutate(class_DMP = case_when(`meth level` == "High" & Direction == "Negative" ~ "A",
                                                                  `meth level` == "Low" & Direction == "Positive" ~ "B",
                                                                  `meth level` == "High" & Direction == "Positive" ~ "C",
                                                                  `meth level` == "Low" & Direction == "Negative" ~ "D",
                                                                  `meth level` == "Intermediate" & `avg MF%` >= 0.5 & Direction == "Positive" ~ "E",
                                                                  `meth level` == "Intermediate" & `avg MF%` <= 0.5 & Direction == "Negative" ~ "F" ,
                                                                  `meth level` == "Intermediate" & `avg MF%` >0.5 & Direction == "Negative" & `avg MF% old` > 0.45 ~ "G",
                                                                  `meth level` == "Intermediate" & `avg MF%` <0.5 &Direction == "Positive" & `avg MF% old` < 0.55 ~"H",
                                                                  `meth level` == "Intermediate" & `avg MF%` <0.55 & Direction == "Negative" & `avg MF% old` < 0.45 ~ "I",
                                                                  `meth level` == "Intermediate" & `avg MF%` >0.45 & Direction == "Positive" & `avg MF% old` > 0.55 ~ "J"))

#create pie chart 
pie <- BIOS_meta_DMPs %>% group_by(`meth level`, Direction) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))

library(ggplot2)
install.packages("webr")
library(webr)
library(dplyr)
pie <- as_tibble(pie) 
names(pie)[names(pie) == "meth level"] <- "MF"

PieDonutCustom <- function (data, mapping, start = getOption("PieDonut.start", 
                                                             0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
                            showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
                                                                                                     0.02), labelposition = getOption("PieDonut.labelposition", 
                                                                                                                                      2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 
                                                                                                                                                                                       0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 
                                                                                                                                                                                                                                              1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                            color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                            showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                            pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                            explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                            family = getOption("PieDonut.family", ""), palette_name="PuBuGn")
{
  (cols = colnames(data))
  if (use.labels) 
    data = moonBook::addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- moonBook::getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = moonBook::getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "x"))
  (donuts = moonBook::getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                             "\n(", scales::percent(df$ratio), ")"), 
                      as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  mainCol = RColorBrewer::brewer.pal(nrow(df), name=palette_name)
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$ratio, 
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                             r0 = as.character(r0), r = as.character(r1), start = "start1", 
                                             end = "end1", fill = pies), alpha = pieAlpha, color = color, 
                                  data = df) + transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df) + geom_text(aes_string(x = "segxend", 
                                                              y = "segyend", label = "label", hjust = "hjust", 
                                                              vjust = "vjust"), size = pieLabelSize, data = df, 
                                                   family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df[df$ratio < labelpositionThreshold, ]) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", 
                           vjust = "vjust"), size = pieLabelSize, 
                data = df[df$ratio < labelpositionThreshold, 
                ], family = family) + geom_text(aes_string(x = "labelx", 
                                                           y = "labely", label = "label"), size = pieLabelSize, 
                                                data = df[df$ratio >= labelpositionThreshold, ], 
                                                family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no", 
                                                 explode = "focus"), alpha = donutAlpha, 
                                      color = color, data = df3)
    }
    else {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no"), 
                                      alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3) + geom_text(aes_string(x = "segxend", 
                                                                 y = "segyend", label = "label", hjust = "hjust", 
                                                                 vjust = "vjust"), size = donutLabelSize, 
                                                      data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", 
                                      y = "labely", label = "label"), size = donutLabelSize, 
                           data = df3, family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3[df3$ratio1 < labelpositionThreshold, 
                              ]) + geom_text(aes_string(x = "segxend", 
                                                        y = "segyend", label = "label", hjust = "hjust", 
                                                        vjust = "vjust"), size = donutLabelSize, 
                                             data = df3[df3$ratio1 < labelpositionThreshold, 
                                             ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                        y = "labely", label = "label"), size = donutLabelSize, 
                                                                             data = df3[df3$ratio1 >= labelpositionThreshold, 
                                                                             ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, 
                          label = title, size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, 
                          y = r3, label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid::grid.newpage()
    print(p1, vp = grid::viewport(height = 1, width = 1))
    print(p3, vp = grid::viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}

setwd("/pvol/EWAS/Blood/Meta-analysis/DMPs")
tiff("Pie chart of DMP classification with direction of change.tiff",
     unit = 'in',
     height = 5,
     width = 5,
     res= 400)
PieDonutCustom(pie, aes(MF,Direction, count=n), 
               #title = "Direction of change with age in DMPs grouped by methylation fraction", 
               color = "white", 
               explode = c(1,2,3), 
               labelposition = 1,
               r0 = 0)
dev.off()

#methylation plot of highly methylated DMP that trends toward 50% with age
cpg <-"cg10501210"
pheno_with_meth <- cbind(BIOS_pheno, meth= as.numeric(BIOS_B[cpg,]))

tiff("High MF entropic.tiff",
     unit = 'in',
     height = 5,
     width = 7,
     res= 400)
ggplot(pheno_with_meth, aes(x=age, y=meth*100)) +
  geom_jitter(width = 0.2, colour = "grey", alpha = 0.95)+
  labs(y=paste("MF%")) +
  geom_smooth(method = "lm", colour = "slategrey") +
  ylim(0,100) + 
  # ggtitle("DMP only") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))
dev.off()


#methylation plot of lowly methylated DMP that trends toward 50% with age

cpg <-"cg04875128"
pheno_with_meth <- cbind(BIOS_pheno, meth= as.numeric(BIOS_B[cpg,]))
tiff("Low MF entropic.tiff",
     unit = 'in',
     height = 5,
     width = 7,
     res= 400)
ggplot(pheno_with_meth, aes(x=age, y=meth*100)) +
  geom_jitter(width = 0.2, colour = "grey", alpha = 0.95)+
  labs(y=paste("MF%")) +
  geom_smooth(method = "lm", colour = "slategrey") +
  ylim(0,100) + 
  # ggtitle("DMP only") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))
dev.off()

#### Cell type adjustment ####
library(tidyverse)

#repeat the formatting of files for the meta-analysis using the cell-type corrected files and format
#according to what had been done for the meta-analysis not adjusted for cell types (i.e. filtering > 3 studied, 
#pvalue adjustment, annotation)

#upload the meta-analysis results for the cell type corrected analysis 
Meta_DMP_CTC <- read_tsv("Meta-analysis results DMPs CTC blood.txt")
#upload the meta-analysis results for the DMP meta-analysis NOT adjusted for cell types 
Meta_DMP <- read_tsv("Meta-analysis results DMPs blood.txt")

#run correlation between the effect size for the cell type adjusted and unadjusted meta-analyses
corr <-cor.test(Meta_DMP$'Effect size', Meta_DMP_CTC$'Effect size', method = "pearson")
#0.94
corr

#Create a new table with only the CpG name, Effect size and P.value for the unadjusted analysis for plotting
ES_plot <- Meta_DMP %>% dplyr::select(CpG, 'Effect size', 'P-value')
#select the same columns in the CTC meta-analysis 
ES_plot_CTC <- Meta_DMP_CTC %>% dplyr::select(CpG, 'Effect size', 'P-value')

#join plots and graph
ES_plot <- ES_plot %>% left_join(ES_plot_CTC, by = "CpG") 

names(ES_plot)
tiff("Corrplot effect size meta-analysis with and without CTC in blood.tiff",
     height = 4,
     width = 7,
     res = 400,
     unit = "in")
ggplot(data = ES_plot)+
  geom_point(aes(x = "Effect size.x", y = "Effect size.y")) +
  ylab("ES meta-analysis adjusted") +
  xlab("ES meta-analysis unadjusted") +
  theme_minimal()+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
dev.off()

