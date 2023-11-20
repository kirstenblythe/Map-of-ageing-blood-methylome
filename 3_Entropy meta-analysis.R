setwd("working directory")

library(metafor)
library(tidyverse)

#### Entropy calculation example ####

#Entropy is calculated in each individual dataset as follows 

#example of how we calculated entropy in each analysis. The JHS entropy calculations are in the 'JHS analysis.R' file provided. calculation 
#load the data for a single dataset e.g. BIOS

B <- data.table::fread("BIOS Normalized beta values after batch correction.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B %>% select(-V1)
M <- logit2(BIOS_B)
pheno <- read.delim("BIOS Phenotypes.txt")

#run limma and adjust for confounders EXCLUDING age on the M values
design = model.matrix(~sex,
                      pheno)
fit1 <- lmFit(M,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M)
plotDensities(M_adj, legend = F)

B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy <- apply(B_adj, 2, calculate_entropy)
Entropy <- as.data.frame(Entropy)
Entropy$Sample_Name <- rownames(Entropy)
pheno <- left_join(pheno, Entropy, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy",
       ylab = "entropy")+
  theme_minimal() 
dev.off()

#run the regression and copy the summary statistics (pval, effect and stderr) to excel file to meta-analyse
fit <- lm(Entropy ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#5.624614e-12 
effect = summary(fit)$coefficients[,1][2]
effect
#effect: 0.0001088927  
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr = 1.567108e-05 

#repeat the entropy calcs using cell type adjusted data
#Entropy with CTC 

designCTC = model.matrix(~sex +
                           CD4T +
                           Bcell +
                           Mono +
                           NK +
                           Gran,
                         pheno)
fit1_CTC <- lmFit(M,
                  designCTC)
fit2_CTC <- eBayes(fit1_CTC)

resid_M_CTC <- residuals(fit2_CTC, M)
M_CTC_adj <- resid_M_CTC + rowMeans(M)
B_CTC_adj <- ilogit2(M_CTC_adj)
plotDensities(B_CTC_adj, legend = F)

Entropy_CTC <- apply(B_CTC_adj, 2, calculate_entropy)
Entropy_CTC <- as.data.frame(Entropy_CTC)
Entropy_CTC$Sample_Name<- rownames(Entropy_CTC)
pheno <- left_join(pheno, Entropy_CTC, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_EntropyCTC.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_CTC)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy with cell type correction",
       ylab = "entropy")+
  theme_minimal()
dev.off()


fitCTC <- lm(Entropy_CTC ~ age, pheno)
summary(fitCTC)
pvalCTC = summary(fitCTC)$coefficients[,4][2]
pvalCTC
#p-value: 1.667965e-15   
effectCTC = summary(fitCTC)$coefficients[,1][2]
effectCTC
#effect: 0.0001098711  
stderrCTC = summary(fitCTC)$coefficients[,2][2]
stderrCTC
#stderr = 1.363851e-05 

cor.test(pheno$Entropy, pheno$Entropy_CTC, method = "pearson")
#pval =   < 2.2e-16
#cor = 0.8784908  

tiff('BIOS Entropy ~ Entropy with cell type correction.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(pheno) +
  geom_jitter(aes(x = Entropy, y = Entropy_CTC), colour = "cornflowerblue") +
  geom_smooth(method = "lm", aes(x = Entropy, y = Entropy_CTC), colour = "black") +
  labs(title = "Pearson's corr = 0.8784908",
       ylabs = "Entropy with cell type correction")+
  theme_minimal()
dev.off()

write.table(pheno, "BIOS Phenotypes.txt",
            sep = "\t",
            col.names = TRUE)

#repeat for ALL age-related CpGS vs non-age-related CpGs

B <- data.table::fread("BIOS Normalized beta values after batch correction.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B %>% select(-V1)

pheno <- read.delim("BIOS Phenotypes.txt", sep = "")

DMPs <- read.csv("/pvol/EWAS/Blood/Meta-analysis/DMPs/Blood_DMPs.csv")
VMPs <- read.csv("/pvol/EWAS/Blood/Meta-analysis/VMPs/Blood_VMPs.csv")
All_cps <- full_join(DMPs, VMPs, by = "MarkerName")
All_cps <- All_cps$MarkerName

B <- as.data.frame(B)
B_all <- B[rownames(B) %in% All_cps,]
M_all <- logit2(B_all)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~predictedSex,
                      pheno)

fit1 <- lmFit(M_all,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M_all)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M_all)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_all <- apply(B_adj, 2, calculate_entropy)
Entropy_all <- as.data.frame(Entropy_all)
Entropy_all$Sample_Name <- rownames(Entropy_all)
pheno$Sample_Name <- as.character(pheno$Sample_Name)
pheno <- left_join(pheno, Entropy_all, by = "Sample_Name")


ggplot(data = pheno, aes(x = age, y = Entropy_all)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy all age-related CpGs",
       ylab = "entropy")+
  theme_minimal()

#calculate entropy on cpgs that do not change with age 

B_none <- B[!rownames(B) %in% All_cps,]
M_none <- logit2(B_none)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~predictedSex,
                      pheno)

fit1 <- lmFit(M_none,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M_none)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M_none)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_none <- apply(B_adj, 2, calculate_entropy)
Entropy_none <- as.data.frame(Entropy_none)
Entropy_none$Sample_Name <- rownames(Entropy_none)
pheno <- left_join(pheno, Entropy_none, by = "Sample_Name")

ggplot(data = pheno, aes(x = age, y = Entropy_none)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy all non-age-related CpGs",
       ylab = "entropy")+
  theme_minimal()

#summary statistics all 
fit <- lm(Entropy_all ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#p-value: 2.484315e-28       
effect = summary(fit)$coefficients[,1][2]
effect
#effect: 0.0002303725      
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr: 2.041515e-05

#summary statistics non-age-related cpgs
fit <- lm(Entropy_none ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#p-value: 0.0005125831       
effect = summary(fit)$coefficients[,1][2]
effect
#effect: -5.696658e-05     
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr: 1.635944e-05 

B <- data.table::fread("BIOS Normalized beta values after batch correction.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B %>% select(-V1)

pheno <- read.delim("BIOS Phenotypes.txt", sep = "")


#### DMPs, DMPs-VMPs, constant VMPs 
DMPs <- read.csv("C:/Users/s4641692/OneDrive - Victoria University/Project/Meta-analysis/Blood/Meta-analysis/DMPs/Blood_DMPs.csv")
VMPs <- read.csv("C:/Users/s4641692/OneDrive - Victoria University/Project/Meta-analysis/Blood/Meta-analysis/VMPs/Blood_VMPs.csv")

DMPs_only <- anti_join(DMPs, VMPs, by = "MarkerName")
DMP_VMP <- inner_join(DMPs, VMPs, by = "MarkerName")
VMPs_only <- anti_join(VMPs, DMPs, by = "MarkerName")

DMPs_only <- DMPs_only$MarkerName
B_DMP <- as.data.frame(B) 
B_DMP <- B_DMP[rownames(B_DMP) %in% DMPs_only,]
M <- logit2(B_DMP)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~predictedSex,
                      pheno)

fit1 <- lmFit(M,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_DMPsonly <- apply(B_adj, 2, calculate_entropy)
Entropy_DMPsonly <- as.data.frame(Entropy_DMPsonly)
Entropy_DMPsonly$Sample_Name <- rownames(Entropy_DMPsonly)
pheno <- left_join(pheno, Entropy_DMPsonly, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_DMPsonly.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_DMPsonly)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy DMPs only",
       ylab = "entropy")+
  theme_minimal() 
dev.off()


fit <- lm(Entropy_DMPsonly ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#0.00263494  
effect = summary(fit)$coefficients[,1][2]
effect
#effect: 4.605342e-05   
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr = 1.528603e-05

#Entropy with CTC 

designCTC = model.matrix(~predictedSex +
                           CD4T +
                           Bcell +
                           Mono +
                           NK +
                           Gran,
                         pheno)
fit1_CTC <- lmFit(M,
                  designCTC)
fit2_CTC <- eBayes(fit1_CTC)

resid_M_CTC <- residuals(fit2_CTC, M)
M_CTC_adj <- resid_M_CTC + rowMeans(M)
B_CTC_adj <- ilogit2(M_CTC_adj)
plotDensities(B_CTC_adj, legend = F)

Entropy_DMPsonly_CTC <- apply(B_CTC_adj, 2, calculate_entropy)
Entropy_DMPsonly_CTC <- as.data.frame(Entropy_DMPsonly_CTC)
Entropy_DMPsonly_CTC$Sample_Name<- rownames(Entropy_DMPsonly_CTC)
pheno <- left_join(pheno, Entropy_DMPsonly_CTC, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_DMPsonly_CTC.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_DMPsonly_CTC)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy DMPs only with cell type correction",
       ylab = "entropy")+
  theme_minimal()
dev.off()


fitCTC <- lm(Entropy_DMPsonly_CTC ~ age, pheno)
summary(fitCTC)
pvalCTC = summary(fitCTC)$coefficients[,4][2]
pvalCTC
#p-value: 0.001221366   
effectCTC = summary(fitCTC)$coefficients[,1][2]
effectCTC
#effect: 4.650944e-05    
stderrCTC = summary(fitCTC)$coefficients[,2][2]
stderrCTC
#stderr =1.43528e-05

cor.test(pheno$Entropy_DMPsonly, pheno$Entropy_DMPsonly_CTC, method = "pearson")
#pval =   < 2.2e-16
#cor = 0.94 


write.table(pheno, "BIOS Phenotypes.txt",
            sep = "\t",
            col.names = TRUE)


VMPs_only <- VMPs_only$MarkerName
B <- as.data.frame(B) 
B_VMP <- B[rownames(B) %in% VMPs_only,]
M <- logit2(B_VMP)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~predictedSex,
                      pheno)

fit1 <- lmFit(M,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_VMPsonly <- apply(B_adj, 2, calculate_entropy)
Entropy_VMPsonly <- as.data.frame(Entropy_VMPsonly)
Entropy_VMPsonly$Sample_Name <- rownames(Entropy_VMPsonly)
pheno <- left_join(pheno, Entropy_VMPsonly, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_VMPsonly.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_VMPsonly)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy VMPsonly",
       ylab = "entropy")+
  theme_minimal() 
dev.off()


fit <- lm(Entropy_VMPsonly ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#1.181521e-12  
effect = summary(fit)$coefficients[,1][2]
effect
#effect: -8.998213e-05    
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr = 1.254433e-05 

#Entropy with CTC 

designCTC = model.matrix(~predictedSex +
                           CD4T +
                           Bcell +
                           Mono +
                           NK +
                           Gran,
                         pheno)
fit1_CTC <- lmFit(M,
                  designCTC)
fit2_CTC <- eBayes(fit1_CTC)

resid_M_CTC <- residuals(fit2_CTC, M)
M_CTC_adj <- resid_M_CTC + rowMeans(M)
B_CTC_adj <- ilogit2(M_CTC_adj)
plotDensities(B_CTC_adj, legend = F)

Entropy_VMPsonly_CTC <- apply(B_CTC_adj, 2, calculate_entropy)
Entropy_VMPsonly_CTC <- as.data.frame(Entropy_VMPsonly_CTC)
Entropy_VMPsonly_CTC$Sample_Name<- rownames(Entropy_VMPsonly_CTC)
pheno <- left_join(pheno, Entropy_VMPsonly_CTC, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_VMPsonly_CTC.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_VMPsonly_CTC)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy VMPs only with cell type correction",
       ylab = "entropy")+
  theme_minimal()
dev.off()


fitCTC <- lm(Entropy_VMPsonly_CTC ~ age, pheno)
summary(fitCTC)
pvalCTC = summary(fitCTC)$coefficients[,4][2]
pvalCTC
#p-value: 3.840255e-11    
effectCTC = summary(fitCTC)$coefficients[,1][2]
effectCTC
#effect: 0.0003075374    
stderrCTC = summary(fitCTC)$coefficients[,2][2]
stderrCTC
#stderr =0.0003075374   

cor.test(pheno$Entropy_VMPsonly, pheno$Entropy_VMPsonly_CTC, method = "pearson")
#pval =   < 2.2e-16
#cor = 0.92

write.table(pheno, "BIOS Phenotypes.txt",
            sep = "\t",
            col.names = TRUE)


#DMPs and VMPs combined

DMP_VMP <- DMP_VMP$MarkerName
B <- as.data.frame(B) 
B_DMPVMP <- B[rownames(B) %in% DMP_VMP,]
M <- logit2(B_DMPVMP)

#run limma and adjust for confounders EXCLUDING age
design = model.matrix(~predictedSex,
                      pheno)

fit1 <- lmFit(M,
              design)
fit2 <- eBayes(fit1)

resid_M <- residuals(fit2, M)

#add the mean beta value for each cpg to the residuals to get a beta matrix 

M_adj <- resid_M + rowMeans(M)
B_adj <- ilogit2(M_adj)
plotDensities(B_adj, legend = F)

#calcuate entropy 
calculate_entropy <- function(B){
  return(sum(B*log2(B) + (1-B)*log2(1-B))/(length(B)*log2(1/2)))
}

Entropy_DMPVMP <- apply(B_adj, 2, calculate_entropy)
Entropy_DMPVMP <- as.data.frame(Entropy_DMPVMP)
Entropy_DMPVMP$Sample_Name <- rownames(Entropy_DMPVMP)
pheno <- left_join(pheno, Entropy_DMPVMP, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_DMPs and VMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_DMPVMP)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy DMP / VMP",
       ylab = "entropy")+
  theme_minimal() 
dev.off()


fit <- lm(Entropy_DMPVMP ~ age, pheno)
summary(fit)
pval = summary(fit)$coefficients[,4][2]
pval
#1.294306e-41  
effect = summary(fit)$coefficients[,1][2]
effect
#effect: 0.000383472    
stderr = summary(fit)$coefficients[,2][2]
stderr
#stderr = 2.745461e-05  

#Entropy with CTC 

designCTC = model.matrix(~predictedSex +
                           CD4T +
                           Bcell +
                           Mono +
                           NK +
                           Gran,
                         pheno)
fit1_CTC <- lmFit(M,
                  designCTC)
fit2_CTC <- eBayes(fit1_CTC)

resid_M_CTC <- residuals(fit2_CTC, M)
M_CTC_adj <- resid_M_CTC + rowMeans(M)
B_CTC_adj <- ilogit2(M_CTC_adj)
plotDensities(B_CTC_adj, legend = F)

Entropy_DMPVMP_CTC <- apply(B_CTC_adj, 2, calculate_entropy)
Entropy_DMPVMP_CTC <- as.data.frame(Entropy_DMPVMP_CTC)
Entropy_DMPVMP_CTC$Sample_Name<- rownames(Entropy_DMPVMP_CTC)
pheno <- left_join(pheno, Entropy_DMPVMP_CTC, by = "Sample_Name")

#Save entropy plot
tiff('BIOS_Entropy_DMPVMP_CTC.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(data = pheno, aes(x = age, y = Entropy_DMPVMP_CTC)) +
  geom_jitter(colour = "deepskyblue4") +
  geom_smooth(method = "lm", colour = "black") +
  labs(title = "BIOS Entropy DMP/VMP with cell type correction",
       ylab = "entropy")+
  theme_minimal()
dev.off()


fitCTC <- lm(Entropy_DMPVMP_CTC ~ age, pheno)
summary(fitCTC)
pvalCTC = summary(fitCTC)$coefficients[,4][2]
pvalCTC
#p-value: 5.297029e-78    
effectCTC = summary(fitCTC)$coefficients[,1][2]
effectCTC
#effect: 0.0003835164   
stderrCTC = summary(fitCTC)$coefficients[,2][2]
stderrCTC
#stderr =1.924813e-05    

cor.test(pheno$Entropy_DMPVMP, pheno$Entropy_DMPVMP_CTC, method = "pearson")
#pval =   < 2.2e-16
#cor = 0.75

write.table(pheno, "BIOS Phenotypes.txt",
            sep = "\t",
            col.names = TRUE)


library(tidyverse)
pheno_entropy <- pheno %>% select(Sample_Name, age, Entropy_DMPsonly, Entropy_VMPsonly, Entropy_DMPVMP)

pheno_entropy <- pheno_entropy %>% pivot_longer(!c(age,Sample_Name),
                                                names_to = "Entropy type",
                                                values_to = "Entropy value")

names(pheno_entropy)[names(pheno_entropy) == "Entropy type"] <- "Condition"

library(viridis)
tiff('BIOS Entropy DMPs, DMP/VMP and VMPs.tiff',
     width =10,
     height = 5,
     units = 'in',
     res = 400)
#pheno_entropy <- pheno_entropy %>% mutate(Class = factor(Class, levels = c("Entropy_VMPs","Entropy_DMPs","Entropy")))
ggplot(pheno_entropy, aes(x = age, y = `Entropy value`))+ 
  geom_jitter(aes(colour = Condition, alpha = 0.8)) +
  geom_smooth(aes(colour = Condition), method = "lm")+ 
  scale_color_manual(values = c("skyblue2","purple", "navy"))+
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal()+
  ggtitle("BIOS Entropy breakdown")
dev.off()

ggsave('BIOS Entropy DMPs, DMP/VMP and VMPs.tiff')


#do the same graph but after cell type correction
pheno_entropy <- pheno %>% select(Sample_Name, age, Entropy_DMPsonly_CTC, Entropy_VMPsonly_CTC, Entropy_DMPVMP_CTC)

pheno_entropy <- pheno_entropy %>% pivot_longer(!c(age,Sample_Name),
                                                names_to = "Entropy type",
                                                values_to = "Entropy value")

names(pheno_entropy)[names(pheno_entropy) == "Entropy type"] <- "Condition"

library(viridis)
tiff('BIOS_Entropy DMPs, DMP/VMP and VMPs in cell type corrected.tiff',
     width =10,
     height = 5,
     units = 'in',
     res = 400)
#pheno_entropy <- pheno_entropy %>% mutate(Class = factor(Class, levels = c("Entropy_VMPs","Entropy_DMPs","Entropy")))
ggplot(pheno_entropy, aes(x = age, y = `Entropy value`))+ 
  geom_jitter(aes(colour = Condition, alpha = 0.8)) +
  geom_smooth(aes(colour = Condition), method = "lm")+ 
  scale_color_manual(values = c("skyblue2","purple", "navy"))+
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal()+
  ggtitle("BIOS with cell type correction")
dev.off()


### Begin analysis from here ####
#### Meta-analysis ####
#load the entropy results

dat <- readxl::read_xlsx("3_Entropy results_Blood_JHS.xlsx", 
                         sheet = "Blood only")

dat$Effect <- as.numeric(dat$Effect)
dat$Stderr <- as.numeric(dat$Stderr)

#fixed effects meta-analysis of entropy and age for genome-wide set of CpGs

fes <- rma(yi=(Effect*10), sei=(Stderr*10), method = "FE", data=dat)
fes

tiff("Entropy forest plot fixed effects blood.tiff",
     height = 10,
     width = 7,
     res = 600,
     unit = 'in')
forest(fes,
       level = 95,
       digits = 4,
       slab = dat$Dataset,
      # header = "Dataset",
       ilab = dat$N,
       ilab.xpos = -0.0075,
       order = -dat$Effect,
       xlab = "Change in entropy (/10 years)",
       header = TRUE,cex = 0.5,
       psize = 0.7,
       xlim=c(-0.012,0.02),
      colout = "#056CF2", 
       col = "blue")
text(-0.0075,58, "N", cex = 0.6, font = 2)
dev.off()

#Entropy with cell type correction

dat$`Effect CTC` <- as.numeric(dat$`Effect CTC`)
dat$`Stderr CTC` <- as.numeric(dat$`Stderr CTC`)
fesCTC <- rma(yi=`Effect CTC`*10, sei=`Stderr CTC`*10, data=dat, method = "FE")
fesCTC

tiff("Entropy forest plot with CTC.tiff",
     height = 10,
     width = 7,
     res = 600,
     unit = 'in')
forest(fesCTC,
       level = 95,
       digits = 4,
       slab = dat$Dataset,
       #header = "Dataset",
       ilab = dat$N,
       ilab.xpos = -0.0075,
       order = -dat$`Effect CTC`,
       xlab = "Change in entropy (/10 years)",
       header = TRUE,cex = 0.5,
       psize = 0.7,
       xlim=c(-0.012,0.015),
       colout = "#05C7F2",
       col = "blue")
text(-0.0075,58, "N", cex = 0.6, font = 2)
dev.off()

cor.test(dat$Effect, dat$`Effect CTC`, method = "pearson")
#pval =  1.611e-05
#cor = 0.5418541   

tiff('Entropy effect size ~ Entropy effect size with cell type correction.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 300)
ggplot(dat) +
  geom_jitter(aes(x = Effect, y = `Effect CTC`), colour = "cornflowerblue") +
  geom_smooth(method = "lm", aes(x = Effect, y = `Effect CTC`), colour = "black") +
 # labs(title = "Pearson's corr = 0.5271369")+
  theme_minimal()
dev.off() 


#Entropy all age-related cpgs vs non age-related cpgs 
fes_all <- rma(yi=`ALL effect`*10, sei=`ALL stderr`*10, data=dat, method = "FE")
fes_all

fes_none <- rma(yi=`NON effect`*10, sei=`NON stderr`*10, data=dat, method = "FE")
fes_none

tiff("Entropy non age-related cpgs forest plot blood.tiff",
     height = 10,
     width = 6,
     res = 600,
     unit = 'in')
forest(fes_none,
       level = 95,
       digits = 4,
       annotate = FALSE,
       slab =dat$Dataset,
       order = -dat$N,
       xlab = "Change in entropy (/10 years)",
       cex.lab = 0.8,
       header = TRUE,
       ilab = dat$N,
       ilab.xpos =-0.012,
       xlim = c(-0.022,0.014),
       cex = 0.6,
       psize = 0.7,
       col = "#747E7E",
       border = "#747E7E",
       colout = "#747E7E")
text(-0.012,58, "N", cex = 0.6, font = 2)
text(0, 59, "Effect size non age-related CpGs",cex = 0.75, font = 2)
dev.off()


tiff("Entropy all age-related cpgs forest plot.tiff",
     height = 10,
     width = 6,
     res = 600,
     unit = 'in')
forest(fes_all,
       level = 95,
       digits = 4,
       annotate = FALSE,
       slab =rep("",length(dat$'Effect')),
       order = -dat$N,
       xlab = "Change in entropy (/10 years)",
      # header = TRUE,
       cex = 0.6,
       psize = 0.7,
       col = "#A3D8C3",
      border = "#A3D8C3",
      colout = "#63C3A1",
      cex.lab = 0.8)
text(0, 59, "Effect size all age-related CpGs",cex = 0.75, font = 2)
dev.off()

# look at change in magnitude of effect between adjusted and CTC adjusted analysis

library(ggplot2)

dat <- dat[1:56,]
dat$Effect <- as.numeric(dat$Effect)
dat$`Effect CTC` <- as.numeric(dat$`Effect CTC`)

tiff("Effect size comparison before and after cell type adjustment.tiff",
     height = 10,
     width = 6,
     res = 600,
     unit = 'in')
ggplot(dat) +
  geom_segment(aes(x = reorder(Dataset, Effect), xend = Dataset, y = Effect*10, yend = `Effect CTC`*10), colour = "grey") +
  geom_point(aes(x = reorder(Dataset, Effect), y = Effect*10), colour = "#056CF2", size = 3) +
  geom_point(aes(x = reorder(Dataset, Effect), y = `Effect CTC`*10), colour = "#05C7F2", size = 3) +
  coord_flip() +
  theme_minimal() +
  ylab("Change in entropy (/10 years)")+
  xlab("Dataset")
dev.off()


#### Entropic vs anti-entropic ####

#### code in the BIOS folder 