##### Perform EWAS of VMPs and age in each dataset ####
#Use GSE197676 as example

#set directory
path = 'pvol/Preprocessing/Blood/GEO/GSE197676'
setwd(path)

library(tidyverse)
library(limma)

resid <- read.table("GSE197676_M_res.txt")
pheno <- read.delim("GSE197676 Phenotypes.txt")

design2 = model.matrix(~ age, pheno)

#Identify non-normal CpGs and remove them
shapirotest <- apply(as.matrix(resid),1,shapiro.test)
pvals <- sapply(shapirotest,"[[",2)
CpGs_to_keep <- names(pvals[pvals>1e-5])
#Remove all sites with raw p-value < 1e-5 (edited) 
resid <- resid[CpGs_to_keep,] 
#353581 cpgs
resid_B <- ilogit2(resid)

#regression of M squared residuals 
sigma2 <- rowSums(resid^2)/ncol(resid) #variance of the residuals
w <- (resid^2) - sigma2
w <- as.matrix(w)

fit1 <- lmFit(w,
              design2)
fit2 <- eBayes(fit1)

#run second regression with beta values to get the effect size and standard error in beta values
sigma2_B <- rowSums(resid_B^2)/ncol(resid_B) #variance of the residuals
w_B <- (resid_B^2) - sigma2_B
w_B <- as.matrix(w_B)

fit1_B <- lmFit(w_B,
                design2)
fit2_B <- eBayes(fit1_B)

coef = "age"
results <- topTable(fit2,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef = coef,
                      number = Inf, 
                      p.value = 1)

results$logFC <- results_B[rownames(results),"logFC"] #extract logFC in beta values
SE <- fit2_B$sigma * fit2_B$stdev.unscaled #extract SE in beta values
results$SE <- SE[rownames(results),coef]

fitted <- fitted(fit2)
fitted_sqrd <- rowSums(fitted^2)
bp <- ncol(resid)*fitted_sqrd/rowSums(w^2)
df <- fit2$rank -1
chisq_pval = sapply(bp,pchisq,df=df,lower.tail = FALSE)
chisq_pval <- as.data.frame(chisq_pval)
results$chisq_pval <- chisq_pval[rownames(results),]
TestStat <- as.data.frame(bp)
names(TestStat)[names(TestStat) == "bp"] <- "TestStat"
results$TestStat <- TestStat[rownames(results), ]

TestStat["cg08319905",] #check that the columns have aligned correctly

#format file for METAL
results_BP <- tibble(CPG = rownames(results),
                     ALLELE1 = rep(1,nrow(results)),
                     ALLELE2 = rep(2,nrow(results)),
                     TESTSTAT = results$TestStat,
                     N = rep(ncol(w),nrow(w)),
                     COEF = results$logFC,
                     SE = results$SE,
                     PVALUE = results$chisq_pval)

#Save p-value histogram
tiff('GSE197676_pvalhist_VMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results_BP,
       mapping = aes(x = PVALUE, binwidth = 20))+
  labs(title="Distribution of raw p-values for age VMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

#save the formatted file for the meta-analysis in METAL
setwd("/pvol/EWAS/Blood/VMPs")
write.table(results_BP,
            file="GSE197676.tbl",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#repeat for remaining datasets 
#### formatting files for METAL ####

#set working directory
setwd("/pvol/EWAS/Blood/VMPs")

#list files in folder
files <- list.files()[grep(".tbl",list.files())]
files 

#multiply the effect size and coefficient for each file because numbers are small

library(tidyverse)
for (f in files)
{
  file <- read_tsv(f) #Read the file
  f <- sub("\\.tbl", "", f) #Obtain the name of the dataset without the ".tbl" at the end
  print(f) #show which dataset we are currently analysing
  
  file$COEF_CORR <- file$COEF*100 #add a column to the file corresponding to the corrected effect size
  file$SE_CORR <- file$SE*100 #add a column to the file corresponding to the corrected standard error
  write.table(file, #save the file
              file=paste0(getwd(),"/",f,"_corrected.tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

#Write code for METAL
#List tbl files
files <- list.files(pattern="corrected.tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME SAMPLESIZE",
                  "MARKER CPG",
                  "WEIGHTLABEL N",
                  "EFFECT EFFECTSIZE_CORR",
                  "STDERR SE_CORR",
                  "SEPARATOR TAB",
                  "PVALUE PVALUE",
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

#set working directory 
setwd("working directory")

list.files()[grep("METAANALYSIS",list.files())]

Meta_VMP <- read.table("2_METAANALYSIS1.TBL",
                       header = TRUE)

##filter on CpGs that are present in at least 15% of samples(N = ~5000)
Meta_VMP <- Meta_VMP %>% filter(Weight >= 5000)

#pvalue adjustment
Meta_VMP <- Meta_VMP %>% mutate(pval_adj = p.adjust(Meta_VMP$P.value, method = "fdr", n = length(Meta_VMP$P.value)))

#add column for significance
Meta_VMP <- Meta_VMP %>% 
  mutate(Sig = ifelse(pval_adj < 0.005, "FDR < 0.005","Not sig"))

#add column for VMP classification
Meta_VMP <- Meta_VMP %>% 
  mutate(coef = ifelse(Zscore >= 0, "Increase", "Decrease"))
Meta_VMP = Meta_VMP %>% mutate(Classification=ifelse(Sig == "Not sig","Not sig",
                                                     ifelse(coef=="Increase","Increase in variance","Decrease in variance")))

#Filter FDR < 0.005
VMPs <- filter(Meta_VMP, pval_adj < 0.005)

write.csv(VMPs, "2_Blood_VMPs.csv")

#### Forest plots ####

#create a forest plot for the VMPs

#Create empty list
L.names <- Meta_VMP$MarkerName
L <- vector("list", length(L.names))
names(L) <- L.names

#create columns in meta-anaysis table for effect size and standard error and these will be empty for the Meta-analysis

Meta_VMP$Effect <- ""
Meta_VMP$Effect <- as.numeric(Meta_VMP$Effect)
Meta_VMP$StdErr <- ""
Meta_VMP$StdErr <- as.numeric(Meta_VMP$StdErr)

#Initiate list
tib <- tibble(CpG = Meta_VMP$MarkerName,
              Dataset = rep("Meta-analysis",nrow(Meta_VMP)),
              ES = Meta_VMP$`Effect`,
              SE = Meta_VMP$`StdErr`,
              PVAL = Meta_VMP$`P.value`,
              FDR = Meta_VMP$pval_adj,
              Type = rep("Meta-analysis",nrow(Meta_VMP)))
tib <- tib %>%
  mutate_if(is.numeric,
            signif,
            digits = 4)

#Initiate the list from meta-analysis
L <- split(tib, seq(nrow(tib))) #splits each row of the tibble into a list component
names(L) <- Meta_VMP$MarkerName

#load the names of the datasets

library(readxl)
list.files()
datasets <- read_excel("C:/Users/s4641692/OneDrive - Victoria University/Project/Meta-analysis/Dataset names for meta-analysis.xlsx")

#Add each study to each component of the list

setwd("/pvol/EWAS/Blood/VMPs")

for (s in datasets$Dataset)
{
  print(s)
  file <- read_tsv(paste0(s,"_corrected.tbl"))
  
#  #Add FDR to the table
  file <- file %>%
    mutate(FDR = p.adjust(PVALUE))
  
  #Add Effect 
  #Create tibble
  summary <- tibble(CpG = file$CPG,
                    Dataset = rep(s,nrow(file)),
                    ES = file$COEF_CORR,
                    SE = file$SE_CORR,
                    PVAL = file$PVALUE,
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
  
#  #Merge the pieces of the two lists that are in common
  L2 <- L[names(subL)]
  L3 <- L[setdiff(names(L),names(subL))]
  listinter <- Map(bind_rows,
                   L2,
                   subL)
  L <- c(listinter,
         L3)
}

setwd("/pvol/EWAS/Blood/Meta-analysis /VMPs")
saveRDS(L,
        file = "List for meta-analysis forest plot blood VMPs.rds")

L <- readRDS("/pvol/EWAS/Blood/Meta-analysis /VMPs/List for meta-analysis forest plot blood VMPs.rds")

setwd("/pvol/EWAS/Blood/Meta-analysis/VMPs")
#forest plot for VMPs 
library(metafor)
cg <- L[["cg21899500"]]

glimpse(datasets)

Datasets <- datasets %>% select(Dataset,N)
cg <- left_join(cg, Datasets, by = "Dataset")
cg <- filter(cg, Dataset != "Meta-analysis")
forest plot of significant VMP
tiff('cg21899500.tiff',
     width =9,
     height = 11,
     units = 'in',
     res = 600)
forest(cg$ES,
       sei= cg$SE,
       slab = cg$Dataset,
       xlab = "% change in DNAm variance per year",
       psize = 0.5,
       header = "Dataset",
       order = desc(cg$N),
       xlim = c(-0.6, 0.5),
       ilab = cg$N,
       ilab.xpos = -0.35,
       cex = 0.7,
       cex.lab = 0.88)
text(-0.35,51,
     "N",
     cex = 0.7)
dev.off()



#### Methylation plots ####
#Visualise most significant cpg across three independent datasets 

#JHS

JHS_B <- data.table::fread("/pvol/Preprocessing/Blood/JHS/JHS beta.txt")
JHS_B <- as.data.frame(JHS_B)
rownames(JHS_B) <- JHS_B$V1
JHS_B <- JHS_B %>% select(-V1)
JHS_pheno <- read.delim("/pvol/Preprocessing/Blood/JHS/JHS Phenotypes.txt", sep = "")

cpg <-"cg21899500"
pheno_with_meth <- cbind(JHS_pheno, meth= as.numeric(JHS_B[cpg,]))
tiff('JHS_cg21899500.tiff',
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

cpg <-"cg21899500"
pheno_with_meth <- cbind(GSE40279_pheno, meth= as.numeric(GSE40279_B[cpg,]))
tiff('GSE40279_cg21899500.tiff',
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
cpg <-"cg21899500"
pheno_with_meth <- cbind(GSE152026_pheno, meth= as.numeric(GSE152026_B[cpg,]))
tiff('GSE152026_cg21899500.tiff',
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

#### Overlap of DMPs and VMPs ####

library(tidyverse)
VMPs <- read_csv("2_Blood_VMPs.csv")
DMPs <- read_csv("1_Blood_DMPs.csv")
L <- list(DMPs$MarkerName,VMPs$MarkerName)
names(L) <- c("DMPs","VMPs")

library(GeneOverlap)

go.object <- newGeneOverlap(L$DMPs,
                            L$VMPs,
                            genome.size = 691302) #total number of cpgs meta-analysed in VMP meta
go.object

#now perform fisher's exact text
go.object <-testGeneOverlap(go.object)
print(go.object)

#### Annotation ####

#load annotation

annotation <- read.delim("1_AnnotationBlood.txt")
#Change name to probeID
Meta_VMP <- dplyr::rename(Meta_VMP,
                          probeID = MarkerName)

meta_anno <- left_join(Meta_VMP, annotation, by= "probeID")

meta_anno <- meta_anno %>% 
  dplyr::select(probeID,
                CpG_chrm,
                CpG_beg,
                CGIposition,
                E062,
                genesUniq_with_enh,
                Weight,
                Zscore,
                P.value,
                Direction,
                pval_adj,
                Sig)

meta_anno <- meta_anno %>% 
  mutate(genesUniq_with_enh = replace_na(genesUniq_with_enh,""))

colnames(meta_anno) = c("CpG",
                         "Chromosome",
                         "Position (hg38)",
                         "CpG island position",
                         "Chromatin state in primary mononuclear cells from peripheral blood",
                         "Annotated gene(s)",
                         "Weight",
                         "Zscore",
                         "Pvalue",
                         "Direction",
                         "FDR",
                         "Significance")

write.table(meta_anno,
            file="Meta-analysis results VMPs blood.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#### Cell type adjustment ####
library(tidyverse)

#repeat the formatting of files for the meta-analysis using the cell-type corrected files and format
#according to what had been done for the meta-analysis not adjusted for cell types (i.e. filtering > 15% samples, 
#pvalue adjustment, annotation)

#upload the meta-analysis results for the cell type corrected analysis 
Meta_VMP_CTC <- read_tsv("2_Meta-analysis results VMPs CTC blood.txt")

#upload the meta-analysis results for the DMP meta-analysis NOT adjusted for cell types 
Meta_VMP <- read_tsv("2_Meta-analysis results VMPs blood.txt")

#Ensure an equal number of probes in both analyses
common_probes <- Meta_VMP %>% inner_join(Meta_VMP_CTC, by = "CpG") %>% pull(CpG)

Meta_VMP <- Meta_VMP %>% filter(CpG %in% common_probes)
Meta_VMP_CTC <- Meta_VMP_CTC %>% filter(CpG %in% common_probes)

#run correlation between the Zscore for the cell type adjusted and unadjusted meta-analyses
corr <-cor.test(Meta_VMP$Zscore, Meta_VMP_CTC$Zscore, method = "pearson")
#0.91
corr

#Create a new table with only the CpG name, Effect size and P.value for the unadjusted analysis for plotting
ES_plot <- Meta_VMP %>% dplyr::select(CpG, Zscore, Pvalue)
#select the same columns in the CTC meta-analysis 
ES_plot_CTC <- Meta_VMP_CTC %>% dplyr::select(CpG, Zscore, Pvalue)

#join plots and graph
ES_plot <- ES_plot %>% left_join(ES_plot_CTC, by = "CpG") 

tiff("Corrplot effect size VMP meta-analysis with and without CTC in blood.tiff",
     height = 4,
     width = 7,
     res = 400,
     unit = "in")
ggplot(data = ES_plot)+
  geom_point(aes(x = Zscore.x, y = Zscore.y)) +
  ylab("Zscore meta-analysis adjusted") +
  xlab("Zscore meta-analysis unadjusted") +
  theme_minimal()+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
dev.off()
