#### preprocessing datasets ####

#### GOLDN ####

#controlled access dataset from the dbGaP

#### WHI ####
#controlled access dataset from the dbGaP

#### GSE55763 ####

#Get the GEO dataset
library(tidyverse)
library(GEOquery)
GSE <- "GSE55763"
setwd(paste0("/pvol/Preprocessing/Blood/GEO/",GSE))

options(timeout = 100000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
download.file(url=paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                         str_sub(GSE,1,5),"nnn/",GSE,"/soft/",
                         GSE,"_family.soft.gz"),
              destfile=paste0(GSE,"_family.soft.gz"))

data <- getGEO(GSE,
               filename = paste0(GSE,"_family.soft.gz"),
               GSEMatrix = FALSE,
               AnnotGPL = FALSE,
               getGPL = FALSE)
GSM_IDs <- data@header$sample_id

library(GEOquery)
pheno <- rbind()
arrayinfo <- c()
for (g in GSM_IDs)
{
  data <- getGEO(g,
                 destdir = getwd(),
                 GSEMatrix = TRUE,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  
  pheno <- rbind(pheno,
                 data@header$characteristics_ch1)
  arrayinfo <- c(arrayinfo,data@header$title)
}

pheno <- list()
for (g in GSM_IDs)
{
  data <- getGEO(g,
                 destdir = getwd(),
                 GSEMatrix = TRUE,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  
  subL <- list(data@header$characteristics_ch1)
  pheno <- c(pheno,
             subL)
}

test <- pheno
colnames <- lapply(test,
                   str_extract,
                   pattern = ".*(?=:)")
test <- lapply(test,
               str_extract,
               pattern = "(?<=:[[:blank:]]).*")
test <- lapply(test,
               as_tibble)
test <- lapply(test,
               t)
test <- map2(test,
             colnames,
             function(first, second) {
               colnames(first) = second
               return(first)
             })

test <- lapply(test,
               as_tibble)
test <- map2(test,
             GSM_IDs,
             function(first, second) {
               first <- first %>% mutate(`GEO accession` = second)
               return(first)
             })
test <- test %>%
  purrr::reduce(full_join)

#Remove soft files that were automatically downloaded
filestoremove <- list.files(pattern=".soft")
sapply(filestoremove,
       file.remove)

#Add info on Sentrix_ID and Sentrix_Position
Sentrix_ID <- str_extract(arrayinfo,
                          pattern = "(?<=,[[:blank:]]).*(?=_)")
Sentrix_Position <- str_extract(arrayinfo,
                                pattern = "(?<=_).*")
test <- test %>%
  mutate(Sentrix_ID = Sentrix_ID,
         Sentrix_Position = Sentrix_Position)

write.csv(test,
          paste(GSE,"phenotypes.csv"))

#Download raw data
options(timeout = 100000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

download.file(url=paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                         str_sub(GSE,1,5),"nnn/",GSE,"/suppl/",GSE,"_unmethylated_methylated_signal_intensities.txt.gz"),
              destfile=paste0(GSE,"_unmethylated_methylated_signal_intensities.txt.gz"))

#Load raw data
raw <- read.delim(paste0(GSE,"_unmethylated_methylated_signal_intensities.txt.gz"),
                  row.names = 1)
detP <- raw %>% select(contains("Pval"))
meth <- raw %>% select(contains(".Methylated.Signal"))
unmeth <- raw %>% select(contains(".Unmethylated"))
tib <- tibble(Sentrix_ID = str_extract(colnames(meth),
                                       pattern = "(?<=X).*(?=_)"),
              Sentrix_Position = str_extract(colnames(meth),
                                             pattern = "(?<=_).*(?=\\.M)"))
test <- read_csv(paste(GSE,"phenotypes.csv")) %>%
  mutate(Sentrix_ID = as.character(Sentrix_ID))
tib <- left_join(tib,
                 test)
colnames(detP) <- tib$`GEO accession`
colnames(meth) <- tib$`GEO accession`
colnames(unmeth) <- tib$`GEO accession`

library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=as.matrix(meth),
                    Unmeth=as.matrix(unmeth),
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
sexpred = getSex(GRset, cutoff = -2)
sexpred <- sexpred %>%
  as_tibble() %>%
  mutate(`GEO accession` = rownames(sexpred)) %>%
  select(`GEO accession`,predictedSex)
test <- left_join(test,
                  sexpred)
which(test$gender!=test$predictedSex) #4 wrong match
write.csv(test,
          paste(GSE,"phenotypes.csv"))

#Load data
library(ChAMP)
rm(methylset)
test <- read_csv(paste(GSE,"phenotypes.csv"))

#Split data in four to avoid R crashing, and then join the tables afterwards
d <- 1:nrow(test)
chunks <- base::split(d,  cut(seq_along(d), 4, labels = FALSE))

for (i in 1:length(chunks))
{
  subtest <- test[chunks[[i]],]
  filter <- champ.filter(beta = getBeta(RSet[,subtest$`GEO accession`]),
                         pd = subtest,
                         Meth = as.matrix(meth[,subtest$`GEO accession`]),
                         UnMeth = as.matrix(unmeth[,subtest$`GEO accession`]),
                         detP = as.matrix(detP[,subtest$`GEO accession`]),
                         arraytype = "450K")
  write.table(filter$beta,
              file = paste0(GSE," beta filtered ",i,".txt"),
              quote=FALSE,
              row.names=TRUE,
              col.names=TRUE,
              sep='\t')
}


#Load all betas and merge
beta1 <- read.table("GSE55763 beta filtered 1.txt")
beta2 <- read.table("GSE55763 beta filtered 2.txt")
beta3 <- read.table("GSE55763 beta filtered 3.txt")
beta4 <- read.table("GSE55763 beta filtered 4.txt")
CpGs <- Reduce(intersect,
               list(rownames(beta1),
                    rownames(beta2),
                    rownames(beta3),
                    rownames(beta4)))
beta <- cbind(beta1[CpGs,],
              beta2[CpGs,],
              beta3[CpGs,],
              beta4[CpGs,])

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
subbeta <- beta
library(tidyverse)
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  subbeta <- subbeta[setdiff(rownames(subbeta),qcprobes),]
}

#Produce quality control graphs to look at the data
subbeta <- as.matrix(subbeta)
champ.QC(beta = subbeta,
         pheno = test$gender,
         dendrogram = FALSE)

write.table(subbeta,
            file="GSE55763 beta filtered.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#Remove temporary beta files
filestoremove <- list.files(pattern="filtered ")
sapply(filestoremove,
       file.remove)

#Normalization of Type I and Type II probes
subbeta <- read.table("GSE55763 beta filtered.txt")
test <- read_csv(paste(GSE,"phenotypes.csv"))
myNorm <- champ.norm(beta=subbeta)
write.table(myNorm,
            file="GSE55763 beta filtered normalised.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#Remove temporary beta files
filestoremove <- list.files(pattern="filtered.txt")
sapply(filestoremove,
       file.remove)

#Data check
champ.QC(beta = myNorm,
         pheno = test$gender,
         dendrogram = FALSE)

library(ChAMP)
library(sva)
library(minfi)
#Split data in two blocks to avoid R crashing at batch correction step
d <- unique(test %>% pull(Sentrix_ID))
chunks <- base::split(d,  cut(seq_along(d), 2, labels = FALSE))

for (i in 1:length(chunks))
{
  subtest <- test %>% filter(Sentrix_ID %in% chunks[[i]])
  subNorm <- myNorm[,subtest$`GEO accession`]
  M <- logit2(subNorm)
  
  myCombat=ComBat(dat=M, #it outputs an M-value matrix adjusted for batch
                  batch=subtest$Sentrix_ID,
                  mod=NULL)
  #Convert back to beta-values after batch correction to run SVD again
  myCombat=ilogit2(myCombat)
  
  write.table(myCombat,
              file = paste0(GSE," beta filtered normalised batch corrected ",i,".txt"),
              quote=FALSE,
              row.names=TRUE,
              col.names=TRUE,
              sep='\t')
}

#Run ComBat to correct batch effects
test <- read_csv(paste(GSE,"phenotypes.csv"))
myNorm <- read.table("GSE55763 beta filtered normalised.txt")


#Remove temporary beta files
filestoremove <- list.files(pattern="normalised.txt")
sapply(filestoremove,
       file.remove)

#Adjust for position on batch
M <- logit2(myCombat)

myCombat=ComBat(dat=M, #it outputs an M-value matrix adjusted for batch
                batch=test$Sentrix_Position,
                mod=NULL)

myCombat=ilogit2(myCombat)

#### preprocessing

#load the data 

pheno = read.delim("GSE55763/GSE55763 Phenotypes.txt", sep = "") 
beta <- data.table::fread("GSE55763/GSE55763_normalized_betas.txt")
beta <- as.data.frame(beta)
beta <- beta %>%
  select(contains("_R"))
rownames(beta) <- beta$ID_REF
beta = beta %>% select(-ID_REF)
keep = rownames(pheno)
beta_filtered <- beta[,keep] 

beta = as.matrix(beta)
#filter using champ.filter
filter = champ.filter(beta = beta,
                      pd = pheno,
                      arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
pheno = filter$pd
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

row.has.na <- apply(beta, 1, function(x){any(is.na(x))})
beta.filtered <- beta[!row.has.na,]
dim(beta.filtered)

#look at quality control graphs 
glimpse(pheno)
pheno$sex <- as.factor(pheno$sex)

champ.QC(beta = beta.filtered,
         pheno = pheno$sex,
         dendrogram = FALSE,
         resultsDir = "./CHAMP_QCimages")

phenofilter$age = as.numeric(phenofilter$age)
ages_range <- ggplot(data = phenofilter) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

ggsave("GSE55763 age distribution.tiff",scale = 1, dpi = "screen")

pheno %>% group_by(sex) %>% summarise(n = n())
pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

write.table(pheno,
            file="GSE55763 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta.filtered,
            file="GSE55763 beta filtered.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE128235 ####
path = "/pvol/Preprocessing/Blood/GEO/GSE128235"
setwd(path)

library(tidyverse)
library(ChAMP)
library(limma)
library(minfi)
library(data.table)
library(GEOquery)

#### formatting
beta = fread("GSE128235_matrix_normalized.txt")
CpGs <- beta$ID_REF
beta = beta %>% select(-ID_REF) %>% as.matrix()
rownames(beta) <- CpGs
detP <- beta %>% as.data.frame() %>% 
  select(contains("Detection"))
drop <- as.character(names(detP))
beta <- as.data.frame(beta)
beta <- beta[,!drop]

beta = beta[,!names(beta) %in% drop]

#upload the phenotypes 

pheno = read.csv("GSE128235 Phenotypes.csv", row.names = 1)

#### preprocessing 
#options(timeout = 100000000000000000000000)
#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
getGEOSuppFiles("GSE128235",
                fetch_files = TRUE,
                filter_regex = ".tar$")
untar("GSE128235/GSE128235_RAW.tar")

#remove the corrupt sample from the phenotype table 
pheno_sub <- pheno %>% filter(geo_accession !="GSM3668181")
write.csv(pheno_sub,"GSE128235 Phenotypes.csv")

#import 
targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = getBeta(RSet),
                      M = getM(RSet),
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")
plotDensities(filter$beta, legend = F)

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

#check for sexes 
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

#Normalisation of type I and type II probes
myNorm <- champ.norm(beta=beta,
                     arraytype="450K")

#check for batch effects
pheno$Slide <- as.factor(pheno$Slide)
pheno$Array <- as.factor(pheno$Array)
pheno$sex <- as.factor(pheno$sex)
pheno$age <- as.numeric(pheno$age)
pheno$diagnosis <- as.numeric(pheno$diagnosis)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            diagnosis),
          resultsDir="./CHAMP_SVDimages/")

#correct for batch 
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            diagnosis),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

#correct for array
M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Array,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            diagnosis),
          resultsDir="./CHAMP_SVDimages/batch_position_corrected/")

champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range 

ggsave("GSE128235 age distribution.tiff")

pheno$age <- as.numeric(pheno$age)
pheno %>% group_by(sex) %>% tally()
pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

percent_male = sum(pheno$sex == "M")/nrow(pheno) * 100

write.table(pheno,
            file="GSE128235 Phenotypes.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

write.table(myCombat,
            file="GSE128235 beta after normalisation and batch correction.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


#### GSE99624 ####
path = "/pvol/Preprocessing/Blood/GEO/GSE99624"
setwd(path)

library(GEOquery)
library(limma)
library(minfi)
library(ChAMP)
library(wateRmelon)
library(tidyverse)
library(doParallel)


#### bulk download of GEO files

files <- c("GSE99624","GSE67705","GSE56105") #these have all been processed on the 450K array

betafiles <- list()
phenofiles <- list()

for (f in files){
  gset = getGEO(f, GSEMatrix = TRUE, getGPL = FALSE)
  if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                    attr(gset, "names")) else idx <- 1
                                    gset = gset[[idx]]
                                    beta = exprs(gset)
                                    betafiles[[f]] = beta
                                    pheno = pData(gset)
                                    phenofiles[[f]] = pheno
}

#### GSE99624 

#extract beta and pheno files from list 
beta = betafiles[["GSE99624"]]
pheno = phenofiles[["GSE99624"]]

#modify the pheno table and save as csv with the raw IDAT files
glimpse(pheno)
pheno <- as_tibble(pheno)
pheno <- pheno %>% dplyr::select(title,
                                 geo_accession,
                                 supplementary_file,
                                 supplementary_file.1,
                                 `age:ch1`,
                                 `disease state:ch1`,
                                 `donor id:ch1`,
                                 `gender:ch1`,
                                 `tissue:ch1`)

names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="disease state:ch1"] <- "disease state"
names(pheno)[names(pheno)=="donor id:ch1"] <- "donor id"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="tissue:ch1"] <- "tissue"

#create columns for Sextrix ID and Sentrix Position 

pheno = within(pheno, supplementary_file<-data.frame(do.call('rbind', strsplit(as.character(supplementary_file), '_', fixed=TRUE))))
glimpse(pheno)

pheno <- cbind(pheno, pheno$supplementary_file[,2])
names(pheno)[names(pheno) == "pheno$supplementary_file[, 2]"] <- "Sentrix_ID"

pheno <- cbind(pheno, pheno$supplementary_file[,3])
names(pheno)[names(pheno) == "pheno$supplementary_file[, 3]"] <- "Sentrix_Position"

pheno <- pheno %>% dplyr::select(-c(3,4))
pheno <- pheno %>% modify_at(c(1,2,4,5,6,7,8), as.factor)
pheno$age <- as.numeric(pheno$age)

rownames(pheno) <- pheno$geo_accession

#save the phenotype file in the IDATs folder
path = paste0(path,"/GSE99624")
setwd(path)

write.csv(pheno,"GSE99624 Phenotypes.csv")

#extract the raw files to check for sexes
options(timeout = 100000000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
getGEOSuppFiles("GSE99624",
                fetch_files = TRUE,
                filter_regex = ".tar$")
path = paste0(path,"/GSE99624")
setwd(path)
untar("GSE99624_RAW.tar")

targets <- read.metharray.sheet(getwd())
RGSet <- read.metharray.exp(targets = targets)
MSet <- preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
predictedSex <- as.matrix(predictedSex)

targets = cbind(targets,predictedSex)
targets$sex[targets$sex == "Female"] <- "F"
targets$sex[targets$sex == "Male"] <- "M"

targets <- targets %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
targets %>% group_by(sexmatch) %>% tally()

#all sexes match 
path = "/volume/Preprocessing/Blood/GEO/GSE99624"
setwd(path)

#run champ filter on the beta matrix to remove SNPs and XY probes 
library(ChAMP)
library(readxl)
sheets <- excel_sheets("/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = impute$beta
for (s in sheets)
{
  qcprobes <- read_excel('/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

filter = champ.filter(beta = as.matrix(beta),
                      pd = targets,
                      arraytype = "450K")
impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

plotDensities(as.matrix(impute$beta), legend = F)
ggsave("beta distribution.tiff")

celltypes = champ.refbase(beta = impute$beta,
                          arraytype = "450K")
celltypes = celltypes[[2]]

beta = impute$beta
pheno = cbind(impute$pd,celltypes)
rownames(pheno) = pheno$geo_accession

path = "/volume/Preprocessing/Blood/GEO/GSE99624"
setwd(path)
write.table(pheno, "GSE99624 Phenotypes.txt")
write.table(beta, "GSE99624 beta.txt")


#### Preprocess

myLoad <- champ.load(directory = paste0(path,"/GSE99624"),
                     filterXY=TRUE, #remove the sex chromosomes
                     filterSNPs=TRUE)

#Produce quality control graphs to look at the data
champ.QC(beta = myLoad$beta,
         pheno = myLoad$pd$sex)

#Normalization of Type I and Type II probes
myNorm <- champ.norm(beta=myLoad$beta,
                     arraytype="450K")

#Produce quality control graphs to look at the data
champ.QC(beta = myLoad$beta,
         pheno = myLoad$pd$sex,
         resultsDir = "./CHAMP_QCimages")

#Convert Slide as a factor so SVD doesn't understand it as a number
myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)

champ.SVD(beta=myNorm,
          pd=myLoad$pd,
          resultsDir="./CHAMP_SVDimages/")

#Adjust for batch effects using the ComBat function using sva package
library(sva)

#Run ComBat to correct batch effects
M <- logit2(myNorm)
batch=as.factor(myLoad$pd$Slide)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=batch,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=myLoad$pd,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

#Save
write.table(myCombat,
            file="GSE99624 Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE115278 ####
path = "/volume/Preprocessing/Blood/GEO/GSE115278"
setwd(path)

library(GEOquery)
library(tidyverse)
library(limma)
library(ChAMP)
library(minfi)
library(reshape)

#### Downloading raw files
#download the raw data files
options(timeout = 1000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000) 

if (!dir.exists("GSE115278")) {
  library(GEOquery)
  getGEOSuppFiles(GEO = "GSE115278",
                  fetch_files = TRUE,
                  filter_regex = "tar$")
}

#untar the IDATS
path = "/volume/Preprocessing/Blood/GEO/GSE115278/GSE115278"
setwd(path)

untar("GSE115278_RAW.tar")

#create the phenotype files
path = "/volume/Preprocessing/Blood/GEO/GSE115278"
setwd(path)

#### 450 preprocessing 
#create the phenotype file for the 450K array samples
list.files()

pheno1 <- read.csv("pheno1.csv")
rownames(pheno1) <- pheno1$X

glimpse(pheno1)

pheno1 <- dplyr::select(pheno1, 
                        title, 
                        geo_accession, 
                        description, 
                        age.ch1, 
                        array.ch1, 
                        Sex.ch1, 
                        study.ch1, 
                        waist.circumference..cm..ch1)

names(pheno1)[names(pheno1) == "age.ch1"] <- "age"
names(pheno1)[names(pheno1) == "Sex.ch1"] <- "sex"
names(pheno1)[names(pheno1) == "study.ch1"] <- "study"
names(pheno1)[names(pheno1) == "waist.circumference..cm..ch1"] <- "waist circumference cm"


pheno450 <- separate(data = pheno1, col = description, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")
names(pheno450)[names(pheno450) == "geo_accession"] <- "Sample_Name"

glimpse(pheno450)
pheno450 <- pheno450 %>% modify_at(c(2, 3, 4, 7, 8), as.factor)
pheno450$age <- as.numeric(pheno450$age)

path = "/volume/Preprocessing/Blood/GEO/GSE115278/450 IDATS"
setwd(path)

write.csv(pheno450, "GSE115278 phenotypes 450.csv")

library(minfi)

#upload the raw files using minfi to first check for sex

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)

filter = champ.filter(beta = getBeta(RSet),
                      M = getM(RSet),
                      pd =  as.matrix(pData(GRset)),
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno450 = filter$pd
pheno450 = as.data.frame(pheno450)
glimpse(pheno450)

#remove the 11 samples with incorrect sex from beta and pheno
pheno450 <- pheno450 %>% modify_at(c(1,3,5,6,7,9,10,14,15), as.factor)
pheno450$age <- as.numeric(pheno450$age)
pheno450$sex <- as.character(pheno450$sex)

pheno450$sex[pheno450$sex == "Male"] <- "M"
pheno450$sex[pheno450$sex == "Female"] <- "F"
glimpse(pheno450)
pheno450$sex <- as.factor(pheno450$sex)

pheno450 <- pheno450 %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
keep <- pheno450 %>% filter(sexmatch == "Yes")
keep = rownames(keep)

dim(beta)
beta = beta[,keep]
dim(beta)

pheno450 <- pheno450 %>% filter(sexmatch == "Yes")

champ.SVD(beta=beta,
          pd=pheno450%>%select(Slide,
                               Array,
                               age,
                               predictedSex,
                               study),
          resultsDir="./CHAMP_SVDimages/")

myNorm <- champ.norm(beta=beta,
                     arraytype="450K")

#Correct for batch
champ.SVD(beta=myNorm,
          pd=pheno450%>%select(Slide,
                               Array,
                               predictedSex,
                               study,
                               age),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)
pheno450$Slide <- as.factor(pheno450$Slide)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno450$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno450%>%select(Slide,
                               Array,
                               predictedSex,
                               study,
                               age),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

limma::plotDensities(myCombat, legend = F)
path = "/volume/Preprocessing/Blood/GEO/GSE115278"
setwd(path)
ggsave("GSE115278 450K beta distribution after preprocessing.tiff")

ages_range <- ggplot(data = pheno450) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

ggsave("GSE115278 450k age distribution.tiff")

pheno450 %>% group_by(predictedSex) %>% tally()

write.table(pheno450, 
            "GSE115278 450K Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = TRUE,
            sep = '\t')

write.table(myCombat,
            file="GSE115278 450K Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### EPIC preprocessing 

#create the phenotype file for the EPIC array samples

path = "/volume/Preprocessing/Blood/GEO/GSE115278"
setwd(path)
list.files()

pheno2 <- read.csv("pheno2.csv")
rownames(pheno2) <- pheno2$X

glimpse(pheno2)

pheno2 <- dplyr::select(pheno2, 
                        title, 
                        geo_accession, 
                        description, 
                        age.ch1, 
                        array.ch1, 
                        Sex.ch1, 
                        study.ch1, 
                        waist.circumference..cm..ch1)

names(pheno2)[names(pheno2) == "age.ch1"] <- "age"
names(pheno2)[names(pheno2) == "Sex.ch1"] <- "sex"
names(pheno2)[names(pheno2) == "study.ch1"] <- "study"
names(pheno2)[names(pheno2) == "waist.circumference..cm..ch1"] <- "waist circumference cm"


phenoEPIC <- separate(data = pheno2, col = description, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")
names(phenoEPIC)[names(phenoEPIC) == "geo_accession"] <- "Sample_Name"

glimpse(phenoEPIC)
phenoEPIC <- phenoEPIC %>% modify_at(c(2, 3, 4, 8), as.factor)
phenoEPIC$age <- as.numeric(phenoEPIC$age)

path = "/volume/Preprocessing/Blood/GEO/GSE115278/EPIC IDATS"
setwd(path)

write.csv(phenoEPIC, "GSE115278 phenotypes EPIC.csv")

library(minfi)


#upload the raw files using minfi to first check for sex

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)

filter = champ.filter(beta = getBeta(RSet),
                      M = getM(RSet),
                      pd =  as.matrix(pData(GRset)),
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

phenoEPIC = filter$pd
phenoEPIC = as.data.frame(phenoEPIC)
glimpse(phenoEPIC)

#remove the samples with incorrect sex from beta and pheno
phenoEPIC <- phenoEPIC %>% modify_at(c(1,3,5,7,8,9,10,15), as.factor)
phenoEPIC$age <- as.numeric(phenoEPIC$age)

glimpse(phenoEPIC)

phenoEPIC$sex[phenoEPIC$sex == "Male"] <- "M"
phenoEPIC$sex[phenoEPIC$sex == "Female"] <- "F"

phenoEPIC$sex <- as.factor(phenoEPIC$sex)

phenoEPIC <- phenoEPIC %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
keep <- phenoEPIC %>% filter(sexmatch == "Yes")
keep = rownames(keep)

dim(beta)
beta = beta[,keep]
dim(beta)

phenoEPIC <- phenoEPIC %>% filter(sexmatch == "Yes")
champ.QC(beta = as.matrix(beta),
         pheno = phenoEPIC$sex)


myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")

#Correct for batch
champ.SVD(beta=myNorm,
          pd=phenoEPIC%>%select(Slide,
                                Array,
                                predictedSex,
                                study,
                                age),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)
phenoEPIC$Slide <- as.factor(phenoEPIC$Slide)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=phenoEPIC$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=phenoEPIC%>%select(Slide,
                                Array,
                                predictedSex,
                                study,
                                age),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

limma::plotDensities(myCombat, legend = F)
path = "/volume/Preprocessing/Blood/GEO/GSE115278"
setwd(path)
ggsave("GSE115278 EPIC beta distribution after preprocessing.tiff")

ages_range <- ggplot(data = phenoEPIC) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE115278 EPIC age distribution.tiff")

phenoEPIC %>% group_by(predictedSex) %>% tally()

write.table(phenoEPIC, 
            "GSE115278 EPIC Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = FALSE,
            sep = '\t')

write.table(myCombat,
            file="GSE115278 EPIC Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE87571 ####

path = "/volume/Preprocessing/Blood/GEO/GSE87571"
setwd(path)

library(GEOquery)
library(tidyverse)
library(limma)
library(minfi)
library(ChAMP)

#### file formatting
#unzip the series matrix file 
series <- read.table("GSE87571_series_matrix.txt", fill = TRUE)
write.csv(series, "Sample characteristics.csv")

#untar the idats
options(timeout = 1000000000000000)
Sys.setenv("VROOM CONNECTION SIZE" = 131072 *10)
getGEOSuppFiles("GSE87571",
                fetch_files = TRUE,
                filter_regex = ".tar$")
path = paste0(path,"/GSE87571")
setwd(path)
untar("GSE87571_RAW.tar")

#read phenotype file and save as csv in IDAT folder 
path = "/volume/Preprocessing/Blood/GEO/GSE87571"
setwd(path)
pheno <- readxl::read_xlsx("GSE87571 Phenotypes.xlsx")

path = paste0(path,"/GSE87571")
setwd(path)
write.csv(pheno, "GSE87571_Phenotypes.csv")

#### Preprocessing

#upload the raw files using minfi to first check for sex
path = "/volume/Preprocessing/Blood/GEO/GSE87571/GSE87571"
setwd(path)

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)

filter = champ.filter(beta = getBeta(RSet),
                      M = getM(RSet),
                      pd =  as.matrix(pData(GRset)),
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")
#Produce quality control graphs to look at the data
limma::plotDensities(filter$beta, legend = F)
ggsave("GSE87571 beta distribution raw.tiff")

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/volume/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = filter$pd
pheno = as.data.frame(pheno)
glimpse(pheno)

#remove the samples with incorrect sex from beta and pheno
pheno <- pheno %>% modify_at(c(1,2,5,18,13,12), as.factor)
pheno$age <- as.numeric(pheno$age)

glimpse(pheno)

pheno$sex <- as.character(pheno$sex)
pheno$predictedSex <- as.character(pheno$predictedSex)

pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()
keep <- pheno %>% filter(sexmatch == "Yes")
keep = rownames(keep)

dim(beta)
beta = beta[,keep]
dim(beta)

pheno <- pheno %>% filter(sexmatch == "Yes")

#Normalization of Type I and Type II probes
myNorm <- champ.norm(beta=beta,
                     arraytype="450K")

#Produce quality control graphs to look at the data after normalization
library(limma)
plotDensities(myNorm, legend = F)
path = "/volume/Preprocessing/Blood/GEO/GSE87571"
setwd(path)
ggsave("GSE87571 beta distribution after normalisation.tiff")

#Singular Value Decomposition (SVD)
#Add phenotypes of interest to identify main sources of variability in DNA methylation overall
#Age, Sex, Batch, Position on the batch, Timepoint, Fitness...
#Convert Slide as a factor so SVD doesn't understand it as a number
pheno$Slide <- as.factor(pheno$Slide)
library(sva)

champ.SVD(beta=myNorm,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/")

#Run ComBat to correct batch effects
M <- logit2(myNorm)
batch=pheno$Slide

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=batch,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

limma::plotDensities(myCombat, legend = F)
path = "/volume/Preprocessing/Blood/GEO/GSE87571"
setwd(path)
ggsave("GSE87571 beta distribution after preprocessing.tiff")

ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range 

ggsave("GSE87571 age distribution.tiff")

Sampleswithage <- subset(pheno, !is.na(age))
keep <- rownames(Sampleswithage)
beta = myCombat[,keep]

ages_range <- ggplot(data = Sampleswithage) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range 
ggsave("GSE87571 age distribution.tiff")
pheno = Sampleswithage

pheno %>% group_by(predictedSex) %>% tally()

write.table(pheno, 
            "GSE87571 Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = FALSE,
            sep = '\t')

write.table(beta,
            file="GSE87571 Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')



#### GSE53740 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE53740"
setwd(path)

library(ChAMP)
library(GEOquery)
library(tidyverse)
library(limma)
library(minfi)

list.files()
pheno <- read.delim("GSE53740 Phenotypes.txt", header = TRUE, sep = " ")
glimpse(pheno)
pheno <- pheno %>% modify_at(c(2:5), as.factor)
pheno <- pheno %>% modify_at(c(1,6,7,8,9,10,11), as.double)

#filter phenotypes on batch 1 

pheno <- pheno %>% filter(batch == 1)

#download the raw files for batch 1 

library(GEOquery)
getGEOSuppFiles("GSE53740",
                fetch_files = TRUE,
                filter_regex = "batch1.csv.gz$")
gunzip("GSE53740/GSE53740_rawSignals_and_detectionPval_batch1.csv.gz")

raw_head <- read.csv("GSE53740/GSE53740_rawSignals_and_detectionPval_batch1.csv", nrow = 5, sep = ,)

list.files()

phenoBatch1 <- read.csv("GSE53740 batch 1 phenotypes with predicted sex.csv")
phenoBatch1 <- phenoBatch1 %>% select(-c(xMed,yMed,predictedSex,sexmatch))
rownames(phenoBatch1) <- phenoBatch1$Sample_Name

BiocManager::install("minfi", lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)
install.packages("rlang",lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)
BiocManager::install("bumphunter",lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)
BiocManager::install("450Kanno.ilmn12.hg19",lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)
BiocManager::install("DMRcate", lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)
install.packages("xfun",lib = "C:/Users/kirst/OneDrive/Documents/R/win-library/4.1",force = TRUE)

library(minfi)

signals <- data.table::fread("GSE53740/GSE53740_rawSignals_and_detectionPval_batch1.csv", sep = ",")
signals <- as.data.frame(signals)
rownames(signals) <- signals$geneInfo
meth <- signals %>% dplyr::select(contains("_MethylatedSignal"))
colnames(meth) <- sub("_MethylatedSignal","",colnames(meth))
meth <- as.matrix(meth)
unmeth <- signals %>% dplyr::select(contains("_UnmethylatedSignal"))
colnames(unmeth) <- sub("_UnmethylatedSignal","",colnames(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% dplyr::select(ends_with("Pval"))
colnames(detP) <- sub("_DetectionPval","",colnames(detP))
detP <- as.matrix(detP)

library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)

beta = getBeta(RSet)
beta <- as.matrix(beta)

library(ChAMP)
filter = champ.filter(beta = beta,
                      pd =  as.matrix(pData(GRset)),
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$Sample_Name <- rownames(pheno)
phenoBatch1$Sample_Name <- sub("X","",phenoBatch1$Sample_Name) 
pheno <- left_join(phenoBatch1, pheno, by = "Sample_Name")

pheno$predictedSex.y <- as.character(pheno$predictedSex.y)

pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex.y,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()
pheno <- pheno %>% filter(sexmatch == "Yes")

pheno %>% summarise(n = n(),
                    range = range(age),
                    sd = sd(age),
                    mean = mean(age))

pheno %>% group_by(sex) %>% tally()

pheno$Sample_Name == colnames(beta)

Sample_Name <- colnames(beta)
Sample_Name <- as.data.frame(Sample_Name)
pheno <- left_join(Sample_Name, pheno)

pheno$Sample_Name == colnames(beta)

write.table(pheno,
            "GSE53740 Phenotypes.txt",
            quote = FALSE,
            row.names = TRUE,
            sep = "\t")

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

pheno <- read.delim("GSE53740 Phenotypes.txt")
pheno <- pheno %>% filter(!is.na(sex))
write.table(pheno,
            "GSE53740 Phenotypes.txt")
beta <- beta[,pheno$Sample_Name]
#Normalization of Type I and Type II probes
beta = as.matrix(beta)
myNorm <- champ.norm(beta=beta,
                     arraytype="450K")


champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

#Singular Value Decomposition (SVD)
#Add phenotypes of interest to identify main sources of variability in DNA methylation overall
#Age, Sex, Batch, Position on the batch, Timepoint, Fitness...
#Convert Slide as a factor so SVD doesn't understand it as a number
pheno$Slide <- as.factor(pheno$Slide)
library(sva)

pheno$Slide <- pheno$Sample_Name
pheno <- separate(pheno, col = Slide, into = c("Slide","Array"), sep = "\\_")

champ.SVD(beta=myNorm,
          pd=pheno %>% select(age,
                              sex,
                              diagnosis,
                              Slide,
                              Array,
                              race,
                              tau_a1.ch1,
                              tau_a2.ch1),
          resultsDir="./CHAMP_SVDimages/")

#Run ComBat to correct batch effects
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno %>% select(age,
                              sex,
                              diagnosis,
                              Slide,
                              Array,
                              race,
                              tau_a1.ch1,
                              tau_a2.ch1),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

myCombat <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno %>% select(age,
                              sex,
                              diagnosis,
                              Slide,
                              Array,
                              race,
                              tau_a1.ch1,
                              tau_a2.ch1),
          resultsDir="./CHAMP_SVDimages/batch_corrected_positioncorrected/")


champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range 

ggsave("GSE53740 age distribution.tiff")

write.table(pheno, 
            "GSE53740 Phenotypes.txt",
            quote = FALSE,
            col.names = TRUE,
            sep = '\t')

write.table(myCombat,
            file="GSE53740 Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### GSE58045 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE58045"
setwd(path)

library(tidyverse)
library(limma)
library(ChAMP)

#### extracting files 
#download beta matrix

gset <- getGEO("GSE58045", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8490", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
plotDensities(ex, legend = F)
write.table(ex,"GSE58045_beta.txt")

beta <- as.matrix(beta)

beta_filtered <- champ.filter(beta = beta, 
                              pd = pheno,
                              filterXY = TRUE,
                              filterSNPs = TRUE,
                              fixOutlier = TRUE)

beta <- beta_filtered[[1]]


beta <- as.matrix(beta)
beta_imputed <- champ.impute(beta = beta,
                             pd = pheno)
which(is.na(beta_imputed))
beta_imputed <- beta_imputed[[1]]

plotDensities(beta_imputed, legend = F)

beta_imputed <- as.matrix(beta_imputed)

write.table(beta_imputed, "GSE58045 beta after filtering.txt")


#extract phenotype table

pheno <- pData(gset)

pheno <- select(pheno,2,38,39,40,43)
glimpse(pheno)

names(pheno)[names(pheno)=="geo_accession"] <- "ID"
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="twin:ch1"] <- "twin"
names(pheno)[names(pheno)=="diseasestate:ch1"] <- "disease"

pheno$ID <- as.factor(pheno$ID)
pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)
pheno$twin <- as.factor(pheno$twin)
pheno$disease <- as.factor(pheno$disease)

#create a column for family numer. This will be used as random intercept (paired design)

pheno <- pheno %>% mutate(family_number = str_extract(twin, "[[:digit:]]+"))
pheno$family_number<- as.factor(pheno$family_number)

write.table(pheno, "GSE58045 Phenotypes.txt")

ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

ggsave("age distribution.tiff")

#### filtering cross reactive probes 
list.files()

beta = read.delim("GSE58045 beta after filtering.txt",
                  header = TRUE,
                  row.names = 1,
                  sep = "")

pheno = read.delim("GSE58045 Phenotypes.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "")

library(readxl)
sheets <- excel_sheets("C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx")
for (s in sheets)
{
  qcprobes <- read_excel('C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno$age = as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE58045 age distribution.tiff",scale = 1, dpi = "screen")

pheno %>% group_by(sex) %>% summarise(n = n())
mean(pheno$age)
sd(pheno$age)

plotDensities(beta, legend = F)
ggsave("GSE77445 beta distribution.tiff",scale = 1, dpi = "screen")

write.table(pheno,
            file="GSE58045 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta,
            file="GSE58045 beta after filtering.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### NAS ####

#controlled access dataset from the dbGaP
#### BIOS ####

#controlled access dataset from the EGA

#### GSE77445 ####
path = "/pvol/Preprocessing/Blood/GEO/GSE77445"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)
library(minfi)

#### downloading datasets
gset <- getGEO("GSE77445",
               GSEMatrix = TRUE,
               getGPL = FALSE)

#the gset will download as a list with one or more elements 

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset2 <- gset[[idx]]

#extract the methylation data matrix 

datamatrix <- exprs(gset2)

#check the number of columns and samples in the data matrix 
dim(datamatrix)

#plot the beta distribution of the data matrix and save 
plotDensities(datamatrix, legend = F)

#save as a beta matrix 
write.table(datamatrix, "GSE77445 beta.txt")

#extract the phenotype file 
pheno <- pData(gset2)

#select only the necessary columns, e.g. ID, age, sex, disease, batch, ethnicity

names(pheno)[names(pheno)=="geo_accession"] <- "ID"
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="bcell proportion:ch1"] <- "Bcell"
names(pheno)[names(pheno)=="granulocytes cell proportion:ch1"] <- "Gran"
names(pheno)[names(pheno)=="monocytes cell proportion:ch1"] <- "Mono"
names(pheno)[names(pheno)=="natural killer cell proportion:ch1"] <- "NK"
names(pheno)[names(pheno)=="cd4 t cell proportion:ch1"] <- "CD4T"
names(pheno)[names(pheno)=="cd8 t cell proportion:ch1"] <- "CD8T"
names(pheno)[names(pheno)=="cortisol stress response area under the curve (auc) with respect to the increase:ch1"] <- "Cortisol stress response AUC"
#save the phenotype file 

write.table(pheno, "GSE77445 Phenotypes.txt")

#### preprocessing 
list.files()

beta = read.delim("GSE77445 beta.txt",
                  header = TRUE,
                  row.names = 1,
                  sep = "")

pheno = read.delim("GSE77445 Phenotypes.txt",
                   header = TRUE,
                   sep = "")

memory.limit(size = 10000000000000000000)
filter = champ.filter(beta = as.matrix(beta),
                      pd = pheno, 
                      arraytype = "450K")

impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = impute$beta
for (s in sheets)
{
  qcprobes <- read_excel('C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = impute$pd
pheno$age = as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE77445 age distribution.tiff",scale = 1, dpi = "screen")

pheno %>% group_by(sex) %>% summarise(n = n())
mean(pheno$age)
sd(pheno$age)

plotDensities(beta, legend = F)
ggsave("GSE77445 beta distribution.tiff",scale = 1, dpi = "screen")

write.table(pheno,
            file="GSE77445 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta,
            file="GSE77445 beta after filtering and imputation.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### JHS ####

#controlled access dataset 

#### FHS ####

#controlled access dataset from the dbGaP

#### GSE49904 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE49904"
setwd(path)


#### downloading
#install GEOquery

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

#install ChAMP

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChAMP")

#install tidyverse

install.packages("tidyverse")

#load GEOquery from the library 

library(GEOquery)
library(ChAMP)
library(tidyverse)

#Download the series matrix file 

gset <- getGEO("GSE49904",
               GSEMatrix = TRUE,
               getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset2 <- gset[[idx]]

#extract the methylation data matrix 

datamatrix <- exprs(gset2)

#check the number of columns and samples in the data matrix 
dim(datamatrix)

#plot the beta distribution of the data matrix and save 
plotDensities(datamatrix, legend = F)

#save as a beta matrix 
write.table(datamatrix, "GSE49904 beta.txt")

#extract the phenotype file 
pheno <- pData(gset2)

#select only the necessary columns, e.g. ID, age, sex, disease, batch, ethnicity
glimpse(pheno)

pheno <- pheno %>% select(1,41,42,43,44,45,46,47,48,49,50)

names(pheno)[names(pheno)=="geo_accession"] <- "Sample_Name"
names(pheno)[names(pheno)=="age (yrs):ch1"] <- "age"
names(pheno)[names(pheno)=="alcohol history:ch1"] <- "alcohol_use"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="height (inch):ch1"] <- "height_inch"
names(pheno)[names(pheno)=="primary diagnosis:ch1"] <- "diagnosis"
names(pheno)[names(pheno)=="tissue:ch1"] <- "tissue"
names(pheno)[names(pheno)=="weight (lb):ch1"] <- "weight_lb"

#save the phenotype file 

write.table(pheno, "GSE100825 Phenotypes.txt")

pheno <- pheno %>% modify_at("age", as.numeric)

ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

#### filtering and imputation
list.files()
beta = read.delim("GSE49904 beta.txt",
                  header = TRUE,
                  row.names = 1,
                  sep = "")

pheno = read.delim("GSE49904 Phenotypes.txt",
                   header = TRUE,
                   sep = "")

memory.limit(size = 10000000000000000000)
library(ChAMP)
filter = champ.filter(beta = as.matrix(beta),
                      pd = pheno)

impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = impute$beta
for (s in sheets)
{
  qcprobes <- read_excel('C:/Users/kirst/OneDrive - Victoria University/DNA methylation/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = impute$pd
pheno$age = as.numeric(pheno$age)

library(tidyverse)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE49904 age distribution.tiff",scale = 1, dpi = "screen")

pheno %>% group_by(sex) %>% summarise(n = n())
mean(pheno$age)
sd(pheno$age)

limma::plotDensities(beta, legend = F)
ggsave("GSE49904 beta distribution.tiff",scale = 1, dpi = "screen")

write.table(pheno,
            file="GSE49904 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta,
            file="GSE49904 beta after filtering and imputation.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#correct for batch

library(sva)
M <- logit2(beta)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$sentrix.barcode.ch1,
                mod=NULL)

myCombat=ilogit2(myCombat)

ChAMP::champ.SVD(beta = myCombat,
                 pd = pheno %>% select(age,
                                       alcohol_use,
                                       sex,
                                       ethnicity.ch1,
                                       diagnosis,
                                       smoking.history.ch1,
                                       tissue,
                                       weight_lb, sentrix.barcode.ch1),
                 resultsDir = "./CHAMP_SVDimages")

write.table(myCombat,
            file="GSE49904 beta after filtering and imputation and batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE80417 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE80417"
setwd(path)

library(GEOquery)
library(tidyverse)
library(data.table)
library(ChAMP)

#### extracting the data files from GEO 
gset <- getGEO("GSE80417", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("13534", attr(gset, "names")) else idx <- 1
gset <- gset[[1]]

pheno_unfiltered <- pData(gset)
write.table(pheno_unfiltered, "GSE80417 Phenotypes raw.txt")

names(pheno)[names(pheno)=="description"] <- "Sample_Name"
row.names(pheno) <- pheno$Sample_Name
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="diseasestatus:ch1"] <- "diseasestatus"
names(pheno)[names(pheno)=="Sex:ch1"] <- "sex"

write.table(pheno, "GSE80417 Phenotypes.txt")

#remove strange people from pheno table 

pheno_filter <- filter(pheno, age != "NA")
pheno_filter <- filter(pheno_filter, age != "891")
pheno_filter <- filter(pheno_filter, age != "883")
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

ggsave("age distribution.tiff")

write.table(pheno_filter, "GSE80417 Phenotypes.txt")

#get the normalised beta values 

GEO <- getGEOSuppFiles("GSE80417", fetch_files = TRUE, filter_regex = "GSE80417_normalizedBetas.csv.gz")

?getGEOSuppFiles

path = paste0(path, "/GSE80417")
setwd(path)
list.files()
gunzip(filename = "GSE80417_normalizedBetas.csv.gz")

#read beta file 
path = paste0(path, "/GSE80417")
setwd(path)
beta <- read_csv("GSE80417_normalizedBetas.csv", col_names = TRUE)
list.files()

rows <- beta$X1

beta <- as.matrix(beta)
row.names(beta) <-beta$X1

beta <- as.data.frame(beta, row.names = TRUE)
dim(beta)

beta[1]
class(beta)
drop <- "X1"
beta = beta[,!(names(beta) %in% drop)]

rownames(beta) <- rows

plotDensities(beta, legend = F)

write.table(beta, "GSE80417 beta.txt")

#download raw files 
options(timeout = 1000000000000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
GEO <- getGEOSuppFiles("GSE80417", fetch_files = TRUE, filter_regex = "rawBetas.csv.gz$")
path = "/pvol/Preprocessing/Blood/GEO/GSE80417/GSE80417"
setwd(path)

untar("GSE80417_RAW.tar")
gunzip("GSE80417_rawBetas.csv.gz")


#check for sex
library(wateRmelon)
beta_raw <- fread("GSE80417_rawBetas.csv", sep = ",")
beta_raw <- as.data.frame(beta_raw)
rownames(beta_raw) <- beta_raw$V1
CpGs <- rownames(beta_raw)
beta_raw <- beta_raw[,-1]
rownames(beta_raw)

path = "/pvol/Preprocessing/Blood/GEO/GSE80417"
setwd(path)
pheno = read.delim("GSE80417 Phenotypes.txt", sep = "")


data("probe.features")
x.probes = filter(probe.features, CHR == "X")
x.probes = rownames(x.probes)

betaX = beta_raw[x.probes,]

library(wateRmelon)
library(matrixStats)
find("colMedians")
betaX <- as.matrix(betaX)
Sex = predictSex(betaX, pc = 2, plot = TRUE)
Sex = as.data.frame(Sex)
Sex$Sample <- rownames(Sex)
names(Sex)[names(Sex) == "Sex"] <- "predictedSex"
pheno$Sample = rownames(pheno)
pheno2 <- merge(pheno, Sex, by = "Sample")
rownames(pheno2) <- pheno2$Sample
pheno2$predictedSex[pheno2$predictedSex == "Male"] <- "M"
pheno2$predictedSex[pheno2$predictedSex == "Female"] <- "F"
pheno2 <- pheno2 %>% mutate(sexmatch = ifelse(sex == predictedSex, "Yes","No"))
pheno2 %>% group_by(sexmatch) %>% tally()
pheno_filter <- pheno %>% filter(sexmatch == "Yes")
#all sexes match so no need to remove bad samples 
path = "/pvol/Preprocessing/Blood/GEO/GSE80417"
setwd(path)
write.table(pheno2, "GSE80417 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

#### filtering and imputation

pheno = read.delim("GSE80417 Phenotypes.txt")
beta = as.matrix(fread("GSE80417 beta.txt"),rownames=1)

#remove the samples without ages from the beta matrix

keep <- rownames(pheno)
beta_filterd <- beta[,keep]

filter = champ.filter(beta = beta_filterd,
                      pd = pheno,
                      arraytype = "450K",
                      filterSNPs = TRUE,
                      filterXY = TRUE)

which(is.na(filter$beta))
champ.QC(beta = filter$beta,
         pheno = filter$pd$sex,
         resultsDir="./CHAMP_QCimages/")

limma::plotDensities(filter$beta)

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno$age = as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE80417 age distribution.tiff",scale = 1, dpi = "screen")

pheno %>% group_by(sex) %>% summarise(n = n())
pheno %>% summarise(mean = mean(age),sd = sd(age))

write.table(pheno,
            file="GSE80417 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta,
            file="GSE80417 beta after filtering and imputation.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### SATSA ####

path = "/pvol/Preprocessing/Blood/ArrayExpress/SATSA"
setwd(path)

#### downloading raw files

#download raw files from ArrayExpress 

library(ArrayExpress)

options(timeout = 1000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*1000)
SATSA = getAE("E-MTAB-7309", type = "raw")

files <- list.files()[grep(".zip$",list.files())]

for (f in files){
  unzip(f)
}

#### preprocessing 

library(ChAMP)
library(GEOquery)
library(tidyverse)
library(limma)
library(minfi)
library(doParallel)
detectCores()
registerDoParallel(cores = 32)

list.files()

path = "/pvol/Preprocessing/Blood/ArrayExpress/SATSA"
setwd(path)

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd = pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")


#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = filter$pd
limma::plotDensities(beta, legend = F)

#Normalization of Type I and Type II probes
myNorm <- champ.norm(beta=beta)

#remove the samples with incorrect sex from beta and pheno
glimpse(pheno)
pheno <- as.data.frame(pheno)

pheno <- pheno %>% 
  modify_at(c(1,18,19), as.factor)
pheno$age <- as.numeric(pheno$age)

pheno$sex[pheno$sex == "male"] <- "M"
pheno$sex[pheno$sex == "female"] <- "F"
glimpse(pheno)

pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
keep <- pheno %>% filter(sexmatch == "Yes")
keep = rownames(keep)

beta = myNorm
beta = beta[,keep]
dim(beta)
pheno <- pheno %>% filter(sexmatch == "Yes")
pheno$twin_pair <- as.factor(pheno$twin_pair)

dim(pheno)

champ.SVD(beta=beta,
          pd=pheno %>%select(age,
                             sex,
                             disease.status,
                             twin.pair,
                             Slide,
                             Array),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(beta)
pheno$Slide <- as.factor(pheno$Slide)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
ChAMP::champ.SVD(beta=myCombat,
                 pd=pheno%>%select(age,
                                   sex,
                                   disease.status,
                                   twin.pair,
                                   Slide,
                                   Array),
                 resultsDir="./CHAMP_SVDimages/batch_corrected/")

M=ilogit2(myCombat)

pheno$Array <- as.factor(pheno$Array)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Array,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
ChAMP::champ.SVD(beta=myCombat,
                 pd=pheno%>%select(age,
                                   sex,
                                   disease.status,
                                   twin.pair,
                                   Slide,
                                   Array),
                 resultsDir="./CHAMP_SVDimages/batch_corrected_position_corrected")

ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("SATSA age distribution.tiff")

pheno %>% group_by(predictedSex) %>% tally()
pheno %>% summarise(mean = mean(age),
                    sd = sd(age))

write.table(pheno, 
            "SATSA Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = TRUE,
            sep = '\t')

write.table(myCombat,
            file="SATSA Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### GSE42861 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE42861"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

options(timeout = 100000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131720 * 100)

gset <- getGEO("GSE42861",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

datamatrix <- exprs(gset)

#check the number of columns and samples in the data matrix 
dim(datamatrix)
plotDensities(datamatrix, legend = F) 

#extract the phenotype files

pheno = pData(gset)
pheno = pheno %>% select(title,
                         geo_accession,
                         supplementary_file,
                         `age:ch1`,
                         `cell type:ch1`,
                         `disease state:ch1`,
                         `gender:ch1`,
                         `smoking status:ch1`,
                         `subject:ch1`)
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="disease state:ch1"] <- "disease state"
names(pheno)[names(pheno)=="subject:ch1"] <- "subject"
pheno <- separate(data = pheno, col = supplementary_file, into = c("Sample_Info","Sentrix_ID", "Sentrix_Position","IDAT colour"), sep = "\\_")

write.csv(pheno,
          "GSE42861 Phenotypes.csv")

#download the raw files 

getGEOSuppFiles("GSE42861",
                fetch_files = TRUE,
                filter_regex = ".tar$")

path = paste0(path,"/GSE42861")
setwd(path)

untar("GSE42861_RAW.tar")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))


filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

#check for sexes

pheno$age <- as.numeric(pheno$age)

pheno$sex[pheno$sex == "m"] <- "M"
pheno$sex[pheno$sex == "f"] <- "F"
glimpse(pheno)
pheno$sex <- as.factor(pheno$sex)
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()


pheno$Slide <- as.factor(pheno$Slide)
pheno$Array <- as.factor(pheno$Array)
pheno$sex <- as.factor(pheno$sex)
pheno$subject <- as.factor(pheno$subject)

champ.SVD(beta=beta,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            subject,
                            smoking.status.ch1),
          resultsDir="./CHAMP_SVDimages/")
myNorm <- champ.norm(beta=beta,
                     arraytype="450K")

#Correct for batch

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            subject,
                            smoking.status.ch1),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

limma::plotDensities(myCombat, legend = F)
path = "/pvol/Preprocessing/Blood/GEO/GSE42861"
setwd(path)
pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE42861 age distribution.tiff")

pheno %>% group_by(predictedSex) %>% tally()
pheno %>% summarise(mean = mean(age),
                    sd = sd(age))

write.table(pheno, 
            "GSE42861 Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = TRUE,
            sep = '\t')

write.table(myCombat,
            file="GSE42861 Normalized beta values after batch correction.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')





#### GSE51032 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE51032"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

options(timeout = 100000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131720 * 100)

gset <- getGEO("GSE51032",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

datamatrix <- exprs(gset)
plotDensities(datamatrix, legend = F) 

pheno = pData(gset)
glimpse(pheno)
pheno <- pheno %>% select(-c(8,9,10,11,12,13,14))
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="cancer type (icd-10):ch1"] <- "cancer type(icd10)"
names(pheno)[names(pheno)=="time to diagnosis:ch1"] <- "time to diagnosis"
pheno <- separate(data = pheno, col = supplementary_file, into = c("Sample_Info","Sentrix_ID", "Sentrix_Position","IDAT colour"), sep = "\\_")

getGEOSuppFiles("GSE51032",
                fetch_files = TRUE,
                filter_regex = ".tar$")
pheno$Sample_Name <- pheno$geo_accession


path = paste0(path,"/GSE51032")
setwd(path)

untar("GSE51032_RAW.tar")

write.csv(pheno,"GSE51032 Phenotypes.csv")

#preprocess using minfi and champ 
?champ.load
import <- champ.load(directory = getwd(),
                     arraytype = "450K")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
#preprocessing failed. Unable to read IDATS? 


#use their preprocessed matrix 

beta = exprs(gset)

#check for sexes 
data("probe.features")
x.probes = probe.features %>% filter(CHR == "X")
x.probes = rownames(x.probes)

betaX = beta[x.probes,]
library(wateRmelon)
predictsex = predictSex(betaX, plot = TRUE)
predictsex <- as.data.frame(predictsex)
predictsex$geo_accession <- rownames(predictsex)
predictsex$predictsex[predictsex$predictsex == "Female"] <- "F"
predictsex$predictsex[predictsex$predictsex == "Male"] <- "M"
pheno <- left_join(pheno,predictsex, by = "geo_accession")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictsex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()
#10 samples do not match the sexes and must be removed. 
rownames(pheno) <- pheno$geo_accession
keep = pheno %>% filter(sexmatch == "Yes")
keep = rownames(keep)
beta = beta[,keep]
pheno = pheno %>% filter(sexmatch == "Yes")

#remove the samples with a cancer diagnosis 

pheno = pheno %>% filter(is.na(`cancer type(icd10)`))
keep <- rownames(pheno)
beta <- beta[,keep]

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      arraytype = "450K")

impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

plotDensities(impute$beta, legend = F)
which(is.na(impute$beta))

beta = impute$beta
pheno = impute$pd

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(pheno)
beta = as.matrix(beta)
glimpse(pheno)
pheno$age <- as.numeric(pheno$age)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE,
         resultsDir="./CHAMP_QCimages")

setwd(path)
pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE51032 age distribution.tiff")

pheno %>% group_by(predictsex) %>% tally()
pheno %>% summarise(mean = mean(age),
                    sd = sd(age))

write.table(pheno, 
            "GSE51032 Phenotypes normal samples only.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = TRUE,
            sep = '\t')

write.table(beta,
            file="GSE51032 filtered and imputed.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')



#### GSE67705 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE67705"
setwd(path)

library(GEOquery)
library(limma)
library(minfi)
library(ChAMP)
library(wateRmelon)
library(tidyverse)
library(doParallel)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*10)
options(timeout = 1000000000000000000000000)
gset = getGEO("GSE67705", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1

gset = gset[[idx]]
pheno = pData(gset)
beta = exprs(gset)

glimpse(pheno)
pheno <- as_tibble(pheno)
pheno <- pheno %>% dplyr::select(title,
                                 geo_accession,
                                 source_name_ch1,
                                 `age:ch1`,
                                 `cell type:ch1`,
                                 `diabetes:ch1`,
                                 `ethnicity:ch1`,
                                 `gender:ch1`,
                                 `tissue:ch1`,
                                 `subject status/id:ch1`)

names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="cell type:ch1"] <- "cell type"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="tissue:ch1"] <- "tissue"
names(pheno)[names(pheno)=="ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno)=="diabetes:ch1"] <- "diabetes"
names(pheno)[names(pheno)=="source_name_ch1"] <- "HIV status"
pheno = as.data.frame(pheno)
rownames(pheno) = pheno$geo_accession

filter = champ.filter(beta = as.matrix(beta),
                      pd = pheno,
                      arraytype = "450K")
impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = impute$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


pheno = impute$pd

pheno$age = as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE66705 age distribution.tiff")

path = "/pvol/Preprocessing/Blood/GEO/GSE67705"
setwd(path)

champ.QC(beta = beta,
         pheno = pheno$`HIV status`,
         dendrogram = FALSE)

champ.SVD(beta=beta,
          pd=pheno %>% select(age,
                              `HIV status`,
                              diabetes,
                              ethnicity),
          resultsDir="./CHAMP_SVDimages/")

write.table(pheno, 
            "GSE67705 Phenotypes.txt",
            quote = FALSE,
            row.names=TRUE,
            col.names = TRUE,
            sep = '\t')

write.table(beta,
            file="GSE67705 beta after filtering and imputation.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### GSE32148 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE32148"
setwd(path)

library(tidyverse)
library(ChaMP)

gset = getGEO("GSE32148", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1

gset = gset[[idx]]
pheno = pData(gset)
beta = exprs(gset)

glimpse(pheno)

pheno <- as_tibble(pheno)
pheno <- pheno %>% dplyr::select(title,
                                 geo_accession,
                                 `age (y):ch1`,
                                 `disease state:ch1`,
                                 `gender:ch1`,
                                 `tissue:ch1`,
                                 `twin:ch1`)
names(pheno_filter)[names(pheno_filter)=="age (y):ch1"] <- "age"
names(pheno)[names(pheno)=="disease state:ch1"] <- "disease state"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"
names(pheno)[names(pheno)=="tissue:ch1"] <- "tissue"
names(pheno)[names(pheno)=="twin:ch1"] <- "twin"

#check sexes using watermelon 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("wateRmelon")

data("probe.features")
x.probes = filter(probe.features, CHR == "X")
x.probes = rownames(x.probes)

betaX = beta[x.probes,]

library(wateRmelon)
Sex = predictSex(betaX, pc = 2, plot = TRUE)
Sex = as.matrix(Sex)

Sex <- cbind(rownames(Sex),Sex)

pheno = cbind(pheno,Sex)

pheno$geo_accession == pheno$`1`
pheno$sex == pheno$`2`

names(pheno)[names(pheno) == "2"] <- "predictedSex"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex, "Yes","No"))

pheno_filter <- pheno %>% filter(sexmatch == "Yes")

keep = rownames(pheno_filter)
beta_filter = beta[,keep]
dim(beta_filter)

plotDensities(beta,legend = F)

which(is.na(beta_filter))

#champ.filter and then champ.impute

filter = champ.filter(beta = as.matrix(beta_filter),
                      pd = pheno_filter,
                      arraytype = "450K")

impute = champ.impute(beta = filter$beta,
                      pd = filter$pd)

plotDensities(impute$beta, legend = F)

which(is.na(impute$beta))

#cell types 

celltypes = champ.refbase(beta = as.matrix(impute$beta),
                          arraytype = "450K")

celltypes = celltypes[[2]]
celltypes = cbind(rownames(celltypes),celltypes)
pheno_filter = cbind(impute$pd,celltypes)
#check that the rows correspond

#check the age distribution

pheno_filter$age = as.numeric(pheno_filter$age)
ages_range <- ggplot(data = pheno_filter) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()

ggsave("GSE32481 age distribution.tiff")

beta = impute$beta
path = "C:/Users/kirst/OneDrive - Victoria University/Data mining/Blood/GSE32148"
setwd(path)

write.table(beta, "GSE32148 beta.txt")
write.table(pheno_filter, "GSE32148 pheno.txt")
getwd()

#filter cross-reactive probes 
beta = read.delim("GSE32148 beta.txt", sep = "")
pheno = read.delim("GSE32148 pheno.txt", sep = "")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno %>% summarise(mean = mean(age),
                    sd = sd(age))

write.table(pheno,
            file="GSE32148 pheno.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(beta,
            file="GSE32148 beta after filtering.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')



#### GSE106648 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE106648"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

options(timeout = 100000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131720 * 100)

gset <- getGEO("GSE106648",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

getGEOSuppFiles("GSE106648",
                fetch_files = TRUE)

path = paste0(path,"/GSE106648")
setwd(path)

gunzip("GSE106648_Matrix_methylated_signal_intensities.txt.gz")
gunzip("GSE106648_Matrix_unmethylated_signal_intensities.txt.gz")

raw <- read.delim("GSE106648_Matrix_methylated_signal_intensities.txt", nrow = 5)

path = paste0(path,"/GSE106648")
setwd(path)

#extract the phenotype file
pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          source_name_ch1,
                          `age:ch1`,
                          `disease status:ch1`,
                          `gender:ch1`,
                          `smoking status:ch1`,
                          `tissue:ch1`)
pheno$preprocessing_info <- pheno$title
pheno <- separate(pheno, col = preprocessing_info, into = c("Sentrix_ID","Sentrix_Position"), sep = "\\_")
names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "disease status:ch1"] <- "disease status"
names(pheno)[names(pheno) == "smoking status:ch1"] <- "smoking status"

meth <- read.delim("GSE106648/GSE106648_Matrix_methylated_signal_intensities.txt", row.names = 1)
unmeth <- read.delim("GSE106648/GSE106648_Matrix_unmethylated_signal_intensities.txt", row.names = 1)

names(meth) <- sub("X","", names(meth))
names(unmeth) <- sub("X","", names(unmeth))
#preprocessing 

#Load annotation to ensure probes on sex chromosomes have been kept
annotation <- read.delim("/pvol/Preprocessing/Annotation.txt")

library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
beta = getBeta(RSet)

#check sex using raw beta 
data(probe.features)

xprobes = probe.features %>% filter(CHR =="X")
xprobes = rownames(xprobes)

betaX = beta[xprobes,]
betaX <- as.matrix(betaX)
predictSex = wateRmelon::predictSex(betaX, plot = TRUE)
predictSex <- as.data.frame(predictSex)
predictSex$title <- rownames(predictSex)
pheno = left_join(pheno, predictSex, by = "title")
pheno$sex[pheno$sex == "male"] <- "Male"
pheno$sex[pheno$sex == "female"] <- "Female"
pheno = pheno %>% mutate(sexmatch = ifelse(sex == predictSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally() #all sexes match 

#Filter probes
library(ChAMP)
filtered <- champ.filter(beta = as.matrix(getBeta(RSet)),
                         M = getM(RSet),
                         pd = pheno,
                         Meth = meth,
                         UnMeth = unmeth,
                         #detP = detP,
                         arraytype = "450K")

champ.QC(beta = filtered$beta,
         pheno = filtered$pd$sex,
         dendrogram = FALSE,
         resultsDir = "./CHAMP_QCimages")

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filtered$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filtered$pd)
glimpse(pheno)
which(is.na(beta))

impute = champ.impute(beta = beta,
                      pd = pheno)

champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE,
         resultsDir = "./CHAMP_QCimages")

#Normalisation of type I and type II probes 
myNorm <- champ.norm(beta = as.matrix(impute$beta))

#Correct for batch 
pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)
pheno$age <- as.numeric(pheno$age)
pheno$`disease status` <- as.factor(pheno$`disease status`)
champ.SVD(beta=myNorm,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            `disease status`,
                            `smoking status`),
          resultsDir="./CHAMP_SVDimages/")

#Run ComBat to correct batch effects
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            `disease status`,
                            `smoking status`),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)
myCombat=ilogit2(myCombat)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            `disease status`,
                            `smoking status`),
          resultsDir="./CHAMP_SVDimages/batch_position_corrected/")

#Produce quality control graphs to look at the data
champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

write.table(myCombat,
            file="GSE106648 beta normalised and batch corrected.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(pheno,
            file="GSE106648 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#### GSE69138 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE69138"
setwd(path)

library(GEOquery)
library(tidyverse)
library(ChAMP)
library(data.table)

#### Downloading datasets 

options(timeout = 100000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

gset <- getGEO("GSE69138",
               GSEMatrix = TRUE,
               getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
pheno_with_ages <- pheno %>% filter(!is.na(`age:ch1`))
getGEOSuppFiles("GSE69138",
                fetch_files = TRUE,
                filter_regex = "185.txt.gz$")

path = "/pvol/Preprocessing/Blood/GEO/GSE69138/GSE69138"
setwd(path)

f <- list.files()

for (f in f){
  gunzip(f)
}

beta <- read.delim("GSE69138_betas_normalized_meneas_185.txt")
meth <- read.delim("GSE69138_meth_raw_meneas_185.txt")
detP <- read.delim("GSE69138_pval_detection_meneas_185.txt")
unmeth <- read.delim("GSE69138_unmeth_raw_meneas_185.txt")
pheno_with_ages$title <- sub("genomic DNA from blood ","X", pheno_with_ages$title)
names(pheno_with_ages)[names(pheno_with_ages) == "title"] <- "Sample_Name"

annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=as.matrix(meth),
                    Unmeth=as.matrix(unmeth),
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRSet <- mapToGenome(RSet)
predictSex <- getSex(GRSet, cutoff = -2)
beta <- getBeta(methylset)
meth <- as.matrix(meth)
unmeth <- as.matrix(unmeth)
detP <- as.matrix(detP)

glimpse(pheno_with_ages)
pheno <- pheno_with_ages %>% select(Sample_Name, 
                                    geo_accession, 
                                    supplementary_file,
                                    `age:ch1`,
                                    `disease state:ch1`,
                                    `gender:ch1`,
                                    `stroke subtype:ch1`)
names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease state"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "stroke subtype:ch1"] <- "stroke subtype"

filter <- champ.filter(beta = beta,
                       pd = pheno,
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP, 
                       arraytype = "450K")

path = "/pvol/Preprocessing/Blood/GEO/GSE69138"
setwd(path)
beta <- filter$beta

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta = beta, arraytype = "450K")

champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

#check sexes
predictSex <- as.data.frame(predictSex)
predictSex$Sample_Name <- rownames(predictSex)

pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"

pheno <- left_join(pheno, predictSex, by = "Sample_Name")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex, "Yes","No")) 
pheno %>% group_by(sexmatch) %>%tally()
pheno <- pheno %>% filter(sexmatch == "Yes")
keep <- pheno$Sample_Name

pheno$age <- as.numeric(pheno$age)
pheno %>% summarise(age = mean(age),
                    sd = sd(age),
                    range = range(age))

beta <- myNorm[,keep]
pheno$Slide <- pheno$Sample_Name
pheno = separate(pheno, col = Slide, into = c("Slide", "Array"), sep = "\\_")
pheno$Slide <- sub("X","", pheno$Slide)

champ.SVD(beta = beta,
          pd = pheno %>% select(age, 
                                sex, 
                                `stroke subtype`,
                                Slide,
                                Array),
          resultsDir = "./CHAMP_SVDimages/")
library(sva)
M <- logit2(beta)
pheno$Slide <- as.factor(pheno$Slide)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(beta = myCombat,
          pd = pheno %>% select(age, 
                                sex, 
                                `stroke subtype`,
                                Slide,
                                Array),
          resultsDir = "./CHAMP_SVDimages/batch_corrected")

write.table(pheno, 
            "GSE69138 Phenotypes.txt",
            col.names = TRUE, 
            sep = "\t")


write.table(myCombat, 
            "GSE69138 beta after normalisation and batch correction.txt",
            col.names = TRUE, 
            row.names = TRUE,
            sep = "\t")

#### GSE50660 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE50660"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

options(timeout = 100000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131720 * 100)

gset <- getGEO("GSE50660",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

getGEOSuppFiles("GSE50660",
                fetch_files = TRUE,
                filter_regex = ".txt.gz$")
?list.files()
files <- list.files(path = paste0(path,"/GSE50660"))
path = "/pvol/Preprocessing/Blood/GEO/GSE50660/GSE50660"
setwd(path)

for (f in files){
  gunzip(f)
}

path = "/pvol/Preprocessing/Blood/GEO/GSE50660"
setwd(path)

beta_processed = read.delim("GSE50660/GSE50660_matrix_processed.txt", row.names = 1, header = TRUE)
pheno = pData(gset)
beta_expressionset <- exprs(gset)
colnames(beta_expressionset)
beta = as.matrix(beta_expressionset)

pheno$Sample_Name <- rownames(pheno)
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
pheno$sex <- as.factor(pheno$sex)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

#check the sexes 
data(probe.features)
xprobes = probe.features %>% filter(CHR =="X")
xprobes = rownames(xprobes)
beta = as.data.frame(beta)
betaX = beta[xprobes,]
betaX <- as.matrix(betaX)
predictSex = wateRmelon::predictSex(betaX, plot = TRUE)
predictSex <- as.data.frame(predictSex)
predictSex$geo_accession <- rownames(predictSex)
pheno = left_join(pheno, predictSex, by = "geo_accession")
pheno = pheno %>% mutate(sexmatch = ifelse(sex == predictSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally() #all sexes match 

#1 sample sex does not match and must be removed from phenotype table and beta matrix 

keep <- pheno %>% filter(sexmatch == "Yes")
rownames(pheno) <- pheno$geo_accession
keep <- rownames(keep)
pheno = pheno[keep,]
beta = beta[,keep]

#Filter probes
library(ChAMP)
filtered <- champ.filter(beta = as.matrix(beta),
                         pd = pheno,
                         arraytype = "450K")


champ.QC(beta = filtered$beta,
         pheno = filtered$pd$sex,
         dendrogram = FALSE)

#remove the cross reactive and SNP probes from Pidsley et al
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filtered$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filtered$pd)

glimpse(pheno)
pheno <- pheno %>% select(title,
                          geo_accession,
                          `age:ch1`,
                          sex,
                          `smoking (0, 1 and 2, which represent never, former and current smokers):ch1`,
                          `tissue:ch1`,
                          Sample_Name,
                          predictSex,
                          sexmatch)
names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "smoking (0, 1 and 2, which represent never, former and current smokers):ch1"] <- "smoking (0, 1 and 2, which represent never, former and current smokers)"
beta = as.matrix(beta)
champ.SVD(beta=beta,
          pd=pheno%>%select(age,
                            sex,
                            `smoking (0, 1 and 2, which represent never, former and current smokers)`),
          resultsDir="./CHAMP_SVDimages/")
write.table(beta,
            file="GSE50660 beta filtered.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(pheno,
            file="GSE50660 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE40279 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE40279"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

list.files()

#load the data 
beta <- read.delim("GSE40279_beta.txt", sep = "")
pheno <- read_csv("GSE40279 pheno with cell types.csv")

#champ.filter
#Replace 0 with smallest positive value
beta[beta == 0] <- min(beta[beta!=0])
#Replace 1 with largest positive value
beta[beta == 1] <- max(beta[beta!=1])

filtered = champ.filter(beta = beta,
                        pd = pheno, 
                        arraytype = "450K",
                        fixOutlier = FALSE)

beta = filtered$beta
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

champ.QC(beta = as.matrix(beta),
         pheno = pheno$sex,
         dendrogram = FALSE)

impute <- champ.impute(beta = as.matrix(beta),
                       pd = pheno)

champ.SVD(beta=impute$beta,
          pd=impute$pd%>%select(age,
                                sex,
                                ethnicity),
          resultsDir="./CHAMP_SVDimages/")

plotDensities(logit2(beta), legend = F)
#there is a small peak where samples with neglibible beta values were included because of no detP value
hist(beta[,3])
quantile(beta[,1], probs = seq(0,1,0.001))[1:10]
x <- beta[,500]
x <- x[x<0.01]
hist(x)

beta_unfiltered<- read.delim("GSE40279_beta.txt", sep = "")
ncol(which(beta_unfiltered <= 0))

beta_unfiltered[beta_unfiltered < 0] <- 0
beta_unfiltered[beta_unfiltered > 1] <- 1

which(beta_unfiltered == 0)

row.is.zero <- apply(beta_unfiltered, 1, function(x){any(x == 0)})
row.is.zero <- row.is.zero[row.is.zero == "FALSE"]
row.is.zero <- as.data.frame(row.is.zero)
keep <- rownames(row.is.zero)

beta_unfiltered <- beta_unfiltered[rownames(beta_unfiltered) %in% keep,]

beta = impute$beta
beta <- beta[rownames(beta) %in% keep,]

plotDensities(logit2(beta), legend = F)

write.table(beta, "GSE40279 beta after filtering and imputation.txt",
            row.names = TRUE,
            col.names = TRUE, 
            sep = "\t")

#### GSE41037 ####

#Get the GEO dataset
directory = "/pvol/Preprocessing/Blood/GEO/GSE41037"
setwd(directory)
library(tidyverse)
library(GEOquery)
data <- getGEO("GSE41037",
               destdir = getwd(),
               GSEMatrix = FALSE,
               AnnotGPL = FALSE,
               getGPL = FALSE)
GSM_IDs <- data@header$sample_id

#Cannot check sex as it is 27K but the MDS plot should tell pretty well if there are two groups that cluster
#Load raw data by downloading data for each sample one by one (the full raw data was not available in ONE big table, but the individual raw data was there on each sample's accession page)
library(GEOquery)
pheno <- rbind()
for (g in GSM_IDs)
{
  data <- getGEO(g,
                 destdir = getwd(),
                 GSEMatrix = FALSE,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  pheno <- rbind(pheno,
                 data@header$characteristics_ch1)
}

test <- pheno
colnames <- as.character(sapply(test[1,],
                                str_extract,
                                pattern = ".*(?=:)"))
colnames(test) <- colnames
test <- test %>%
  as_tibble%>%
  mutate_all(str_extract,
             pattern = "(?<=:[[:blank:]]).*")%>%
  mutate(`GEO accession` = GSM_IDs)%>%
  dplyr::rename(Sentrix_ID = sentrixbarcode,
                Sentrix_Position = sentrixposition)

write_csv(test,
          file = "GSE41037 phenotypes with bad samples.csv")

#Load raw data
pheno <- read_csv("GSE41037 phenotypes with bad samples.csv")
raw <- read.table('GSE41037_non_normalized.txt.gz',
                  header = TRUE)
CpGs <- raw$TargetID
raw <- as.matrix(raw[,-c(1:4)])
rownames(raw) <- CpGs

unmeth <- raw[,grepl(".Signal_A",colnames(raw))]
colnames(unmeth) <- pheno$`GEO accession`

meth <- raw[,grepl(".Signal_B",colnames(raw))]
colnames(meth) <- pheno$`GEO accession`

detP <- raw[,grepl(".Detec",colnames(raw))]
colnames(detP) <- pheno$`GEO accession`


#Now put data into an methylset object
library(minfi)
annotation <- c("IlluminaHumanMethylation27k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)

#Remove samples with > 10% of probes with detection p-value > 0.01
filter=function(v)
{
  thresh=0.1*length(v)
  if (length(which(v>0.01))<thresh)
    return("OK")
  else
    return("Bad sample")
}
badsamples <- which(apply(detP,2,filter)=="Bad sample") #12 samples failed QC
pheno <- pheno[-badsamples,]
RSet <- RSet[,-badsamples]

#Probe filtering
filter2=function(v)
{
  thresh=floor(0.05*nrow(pheno)) #5% of samples
  if (length(which(v>0.01))<thresh)
    return("OK")
  else
    return("Bad probe")
}
p=apply(detP,1,filter2)
length(which(p=="Bad probe"))  #2309 probes did not pass QC
CpGs <- rownames(RSet)
CpGs_tokeep=setdiff(CpGs,CpGs[which(p=="Bad probe")])

#Remove probes with general SNP list
library(ChAMP)
data(EPIC.manifest.pop.hg19)
maskname <- rownames(EPIC.manifest.pop.hg19)[which(EPIC.manifest.pop.hg19$MASK_general_EUR == 
                                                     TRUE)]
CpGs_SNPs=intersect(CpGs_tokeep,maskname)
CpGs_tokeep=setdiff(CpGs_tokeep,CpGs_SNPs)

#Remove list of cross-hybridising probes and SNP probes from Peters et al.
library(readxl)
sheets <- excel_sheets("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data/Cross reactive and SNP probes from Pidsley et al.xlsx")
library(tidyverse)
for (s in sheets)
{
  qcprobes <- read_excel('F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  CpGs_tokeep=setdiff(CpGs_tokeep,qcprobes)
}

#Get beta
B <- getBeta(RSet)

#Fix outliers
B=B[CpGs_tokeep,]
#Replace 0 with smallest positive value
B[B == 0] <- min(B[B!=0])
#Replace 1 with largest positive value
B[B == 1] <- max(B[B!=1])

#Produce quality control graphs to look at the data
champ.QC(beta = B,
         pheno = pheno$gender,
         dendrogram = FALSE)
pheno <- pheno %>%
  mutate(gender = ifelse(is.na(gender),
                         "unknown",
                         gender))
QC.GUI(beta = B,
       pheno = pheno$gender) #There are 5 samples that do not cluster with their respective sex or tissue: GSM1007609, GSM1007668, GSM1007697, GSM1007220, GSM1007207 --> remove them
pheno <- pheno %>%
  dplyr::filter(!`GEO accession`%in%c("GSM1007609", "GSM1007668", "GSM1007697", "GSM1007220", "GSM1007207"))
B <- B[,pheno$`GEO accession`]
#Sample GSM1007327 has no assigned sex but clusters with females --> pretty safe to say she is female
pheno <- pheno %>%
  mutate(gender = ifelse(`GEO accession`=="GSM1007327",
                         "female",
                         gender))

write.table(pheno,
            file="GSE41037_phenotypes2.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Remove sex chr probes
#Load annotation
annotation <- read.delim("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data/Annotation.txt")
CpGs_tokeep <- setdiff(CpGs_tokeep,
                       annotation%>%dplyr::filter(CpG_chrm%in%c("chrX","chrY"))%>%pull(probeID))
B <- B[CpGs_tokeep,]

#Correct for batch
pheno <- pheno %>%
  mutate_at(vars(gender,diseasestatus:`GEO accession`),
            as.factor)
champ.SVD(beta=B,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/")

#Run ComBat to correct batch effects
library(sva)
M <- logit2(B)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

#Run ComBat to correct position on batch
library(sva)
M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_position_corrected/")

write.table(myCombat,
            file="GSE41037_filtered_batch_corrected_beta.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

#Produce quality control graphs to look at the data
champ.QC(beta = myCombat,
         pheno = pheno$gender,
         dendrogram = FALSE)
QC.GUI(beta = B,
       pheno = pheno$gender)


#### GSE41169 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE41169"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

options(timeout = 100000000000000000000)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131720 * 100)

gset <- getGEO("GSE41169",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
beta <- exprs(gset)
phenoraw <-pData(gset)
glimpse(phenoraw)
pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             `age:ch1`,
                             `diseasestatus (1=control, 2=scz patient):ch1`,
                             `gender:ch1`,
                             `plate:ch1`,
                             `sample type:ch1`,
                             `sentrix barcode:ch1`,
                             `sentrix position:ch1`,
                             `well id:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "diseasestatus (1=control, 2=scz patient):ch1"] <- "diseasestatus (1=control, 2=scz)"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "`plate:ch1"] <- "plate"

#check sexes

colnames(beta)
data(probe.features)
x.probes = filter(probe.features, CHR == "X")
x.probes = rownames(x.probes)
predictSex <- wateRmelon::predictSex(beta, x.probes = x.probes, pc = 2, plot = TRUE)
predictSex <- as.data.frame(predictSex)
predictSex$geo_accession <- rownames(predictSex)
pheno <- left_join(pheno, predictSex, geo_accession, by = "geo_accession")
pheno <- mutate(pheno, sexmatch = ifelse(sex == predictSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()
pheno <- pheno %>% filter(sexmatch == "Yes")

keep <- pheno$geo_accession
beta <- beta[,keep]

names(pheno)[names(pheno) == "geo_accession"] <- "Sample_Name"

filter <- champ.filter(beta = as.matrix(beta),
                       pd = pheno, 
                       arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


pheno = as.data.frame(filter$pd)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

impute <- champ.impute(beta = as.matrix(beta),
                       pd = pheno)
champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE)

beta = impute$beta
pheno = impute$pd

champ.SVD(beta=beta,
          pd=pheno%>%select(age,
                            `diseasestatus (1=control, 2=scz)`,
                            sex,
                            `sentrix barcode:ch1`,
                            `sentrix position:ch1`,
                            `well id:ch1`,
                            `plate:ch1`),
          resultsDir="./CHAMP_SVDimages/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE41169 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

write.table(pheno,
            "GSE41169 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(beta,
            "GSE41169 beta after filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE53840 ####

#Get the GEO dataset
library(tidyverse)
directory = "/pvol/Preprocessing/Blood/GEO/GSE53840"
setwd(directory)
library(GEOquery)
data <- getGEO("GSE53840",
               destdir = getwd(),
               GSEMatrix = FALSE,
               AnnotGPL = FALSE,
               getGPL = FALSE)
GSM_IDs <- data@header$sample_id

library(GEOquery)
pheno <- rbind()
for (g in GSM_IDs)
{
  data <- getGEO(g,
                 destdir = getwd(),
                 GSEMatrix = FALSE,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  pheno <- rbind(pheno,
                 data@header$characteristics_ch1)
}

test <- pheno
colnames <- as.character(sapply(test[1,],
                                str_extract,
                                pattern = ".*(?=:)"))
colnames(test) <- colnames
test <- test %>%
  as_tibble%>%
  mutate_all(str_extract,
             pattern = "(?<=:[[:blank:]]).*")%>%
  mutate(`GEO accession` = GSM_IDs)%>%
  dplyr::rename(Sentrix_ID = samplechips12,
                Sentrix_Position = stripe)

write_csv(test,
          file = "GSE53840 phenotypes.csv")

#Check sex
pheno <- read_csv("GSE53840 phenotypes.csv")
memory.limit(1000000000)
setwd(directory)

#Load raw data
raw <- read.table('GSE53840_non_normalized.txt.gz',
                  header = TRUE)
CpGs <- raw$ID_REF
raw <- as.matrix(raw[,-1])
rownames(raw) <- CpGs

unmeth <- raw[,grepl(".Unmeth",colnames(raw))]
colnames(unmeth) <- pheno$`GEO accession`

meth <- raw[,grepl(".Meth",colnames(raw))]
colnames(meth) <- pheno$`GEO accession`

detP <- raw[,grepl(".Detec",colnames(raw))]
colnames(detP) <- pheno$`GEO accession`

#Now put data into an methylset object to check Sex with minfi
library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset,
                       cutoff = -2)$predictedSex #all male

#Filter probes
library(ChAMP)
filtered <- champ.filter(beta = getBeta(RSet),
                         M = getM(RSet),
                         pd = pheno,
                         Meth = meth,
                         filterXY = FALSE,
                         UnMeth = unmeth,
                         detP = detP)

library(readxl)
sheets <- excel_sheets("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data/Cross reactive and SNP probes from Pidsley et al.xlsx")
subbeta <- filtered$beta
library(tidyverse)
for (s in sheets)
{
  qcprobes <- read_excel('F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  subbeta <- subbeta[setdiff(rownames(subbeta),qcprobes),]
}
filtered$beta <- subbeta

#Produce quality control graphs to look at the data
champ.QC(beta = filtered$beta,
         pheno = pheno$hivstatus,
         dendrogram = FALSE)

#Normalization of Type I and Type II probes
myNorm <- champ.norm(beta=filtered$beta)

#Have a look at variables
setwd(directory)
pheno <- pheno%>%
  mutate_at(vars(hivstatus,
                 Sentrix_ID:`GEO accession`),
            as.factor)

library(ChAMP)
champ.SVD(beta=myNorm,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/")

library(sva)

#Run ComBat to correct batch effects
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

#Convert back to beta-values after batch correction to run SVD again
myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)
myCombat=ilogit2(myCombat)

write.table(myCombat,
            file="GSE53840_normalized_batch_corrected_beta.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

champ.SVD(beta=myCombat,
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_position_corrected/")

#Produce quality control graphs to look at the data
champ.QC(beta = myCombat,
         pheno = pheno$hivstatus,
         dendrogram = FALSE)


#### GSE67751 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE67751"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE67751",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
beta <- exprs(gset)
phenoraw <-pData(gset)
glimpse(phenoraw)
pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             `age:ch1`,
                             `ageaccelerationresidual:ch1`,
                             `ageaccelerationvscontrols:ch1`,
                             `b cells:ch1`,
                             `cd4 t cells:ch1`,
                             `cd8 t cells:ch1`,
                             `dnamage:ch1`,
                             `exhausted cd8t cells:ch1`,
                             `gender:ch1`,
                             `hivstatus:ch1`,
                             `monocytes:ch1`,
                             `naïve cd4 t cells:ch1`,
                             `naïve cd8 t cells:ch1`,
                             `natural killer:ch1`,
                             `tissue:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "b cells:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4 t cells:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8 t cells:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "exhausted cd8t cells:ch1"] <- "Exhausted CD8T"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "hivstatus:ch1"] <- "HIVstatus"
names(pheno)[names(pheno) == "monocytes:ch1"] <- "Mono"
names(pheno)[names(pheno) == "naïve cd4 t cells:ch1"] <- "naïve CD4T"
names(pheno)[names(pheno) == "naïve cd8 t cells:ch1"] <- "naïve CD8T"
names(pheno)[names(pheno) == "natural killer:ch1"] <- "NK"

#check sexes
colnames(beta)
data(probe.features)
x.probes = filter(probe.features, CHR == "X")
x.probes = rownames(x.probes)
predictSex <- wateRmelon::predictSex(beta, x.probes = x.probes, pc = 2, plot = TRUE)
predictSex <- as.data.frame(predictSex)
predictSex$geo_accession <- rownames(predictSex)
pheno <- left_join(pheno, predictSex, geo_accession, by = "geo_accession")
pheno$sex[pheno$sex == "female"] <- "Female"
pheno$sex[pheno$sex == "male"] <- "Male"
pheno <- mutate(pheno, sexmatch = ifelse(sex == predictSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

filter <- champ.filter(beta = as.matrix(beta),
                       pd = pheno, 
                       arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


pheno = as.data.frame(filter$pd)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

impute <- champ.impute(beta = as.matrix(beta),
                       pd = pheno)
champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE)

beta = impute$beta
pheno = impute$pd

champ.SVD(beta=beta,
          pd=pheno%>%select(age,
                            HIVstatus,
                            sex),
          resultsDir="./CHAMP_SVDimages/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE67751 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

write.table(pheno,
            "GSE67751 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(beta,
            "GSE67751 beta after filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE72775 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE72775"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE72775",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
phenoraw <-pData(gset)
glimpse(phenoraw)

pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             `age:ch1`,
                             `b cell:ch1`,
                             `bioage4:ch1`,
                             `cd4.naive:ch1`,
                             `cd8.naive:ch1`,
                             `cd4t:ch1`,
                             `cd8.naive:ch1`,
                             `cd8pcd28ncd45ran:ch1`,
                             `cd8t:ch1`,
                             `dnamage:ch1`,
                             `ethnicity:ch1`,
                             `granulocyte:ch1`,
                             `monocyte:ch1`,
                             `natural killer cell:ch1`,
                             `plasmablast:ch1`,
                             `Sex:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "b cell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno) == "granulocyte:ch1"] <- "Gran"
names(pheno)[names(pheno) == "monocyte:ch1"] <- "Mono"
names(pheno)[names(pheno) == "plasmablast:ch1"] <- "Plasmablast"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "natural killer cell:ch1"] <- "NK"
pheno$title <- sub("PEGforRace: genomic DNA from whole blood of subject ","", pheno$title)

getGEOSuppFiles("GSE72775",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE72775/GSE72775"
setwd(path)

gunzip("GSE72775_datSignal.csv.gz")
gunzip("GSE72775_datBetaNormalized.csv.gz")

beta_processed <- read_csv("GSE72775_datBetaNormalized.csv")
signals <- read_csv("GSE72775_datSignal.csv")

#preprocess using signal intensities

rownames(signals) <- signals$TargetID
signals <- as.data.frame(signals)
meth <- signals %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Unmethylated","",colnames(meth))
meth <- meth %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Methylated","", colnames(meth))
unmeth <- signals %>% select(contains("Unmethylated"))
colnames(unmeth) <- sub("Unmethylated","", colnames(unmeth))
detP <- signals %>% select(contains("DectionPvalue"))
colnames(detP) <- sub("DectionPvalue","", colnames(detP))
meth <- as.matrix(meth)
unmeth <- as.matrix(unmeth)
detP <- as.matrix(detP)

annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRSet <- mapToGenome((RSet))
predictSex <- getSex(GRSet, cutoff = -2)
predictSex$title <- rownames(predictSex)
predictSex <- as.data.frame(predictSex)

pheno <- left_join(pheno, predictSex, by = "title")
pheno$sex[pheno$sex == "male"] <- "M"
pheno$sex[pheno$sex == "female"] <- "F"
pheno <- mutate(pheno, sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally() #all sexes match 
pheno$Sample_Name <- pheno$title

filter <- champ.filter(beta = getBeta(RSet),
                       pd = pheno,
                       Meth = meth,
                       UnMeth = unmeth,
                       detP = detP,
                       arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

which(is.na(beta))

#normalisation 

myNorm <- champ.norm(beta = beta)

champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta = myNorm,
          pd = pheno %>% select(age,
                                ethnicity,
                                sex),
          resultsDir = "./ChAMP_SVDimages")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE72775 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

path = "/pvol/Preprocessing/Blood/GEO/GSE72775"
setwd(path)

write.table(pheno,
            "GSE72775 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myNorm,
            "GSE72775 beta after normalisation filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")      

#### GSE11169 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE111629"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE111629",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)
pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          supplementary_file,
                          `age:ch1`,
                          `disease state:ch1`,
                          `gender:ch1`,
                          `ethnicity:ch1`,
                          `tissue:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease state"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "source_name_ch1"] <- "Sample_Name"
pheno$Sample_Name <- sub("X","",pheno$Sample_Name) 
pheno$Sentrix_ID <- pheno$Sample_Name
pheno <- separate(pheno, col = Sentrix_ID, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")

getGEOSuppFiles("GSE111629",
                fetch_files = TRUE,
                filter_regex = "tar$")

untar("GSE111629/GSE111629_RAW.tar")
getwd()

pheno <- pheno %>% filter(supplementary_file != "NONE")
write.csv(pheno,
          "GSE111629 Phenotypes.csv")

targets <- read.metharray.sheet(getwd())
targets <- filter(targets, Slide != "3999984061")
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))


filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="450K")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state,
                            ethnicity),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state,
                            ethnicity),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")


M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state,
                            ethnicity),
          resultsDir="./CHAMP_SVDimages/batch_Position_corrected/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE111629 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()
head(pheno)[,1:5]
write.table(pheno,
            "GSE111629 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE111629 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE72774 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE72774"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE72774",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
phenoraw <-pData(gset)
glimpse(phenoraw)

pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             `age:ch1`,
                             `bcell:ch1`,
                             `bioage4:ch1`,
                             `cd4.naive:ch1`,
                             `cd8.naive:ch1`,
                             `cd4t:ch1`,
                             `cd8.naive:ch1`,
                             `cd8pcd28ncd45ran:ch1`,
                             `cd8t:ch1`,
                             `diseasestatus:ch1`,
                             `dnamage:ch1`,
                             `ethnicity:ch1`,
                             `familyhistoryofparkinsons:ch1`,
                             `gran:ch1`,
                             `levodopaperdayincludingagonists:ch1`,
                             `levodopastatus:ch1`,
                             `mono:ch1`,
                             `nk:ch1`,
                             `numberofyearsinschool:ch1`,
                             `packyearstotal:ch1`,
                             `plasmablast:ch1`,
                             `Sex:ch1`,
                             `smokingstatus:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "b cell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "diseasestatus:ch1"] <- "disease.status"
names(pheno)[names(pheno) == "familyhistoryofparkinsons:ch1"] <- "familyhistoryofparkinsons"
names(pheno)[names(pheno) == "ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno) == "gran:ch1"] <- "Gran"
names(pheno)[names(pheno) == "mono:ch1"] <- "Mono"
names(pheno)[names(pheno) == "plasmablast:ch1"] <- "Plasmablast"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "natural killer cell:ch1"] <- "NK"
pheno$title <- sub("PEGcaucasian: genomic DNA from whole blood of subject ","", pheno$title)

getGEOSuppFiles("GSE72774",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE72774/GSE72774"
setwd(path)

gunzip("GSE72774_datSignal.csv.gz")
gunzip("GSE72774_datBetaNormalized.csv.gz")

beta_processed <- read_csv("GSE72774_datBetaNormalized.csv")
signals <- read_csv("GSE72774_datSignal.csv")

#preprocess using signal intensities

rownames(signals) <- signals$TargetID
signals <- as.data.frame(signals)
meth <- signals %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Unmethylated","",colnames(meth))
meth <- meth %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Methylated","", colnames(meth))
unmeth <- signals %>% select(contains("Unmethylated"))
colnames(unmeth) <- sub("Unmethylated","", colnames(unmeth))
detP <- signals %>% select(contains("DectionPvalue"))
colnames(detP) <- sub("DectionPvalue","", colnames(detP))
meth <- as.matrix(meth)
unmeth <- as.matrix(unmeth)
detP <- as.matrix(detP)

annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRSet <- mapToGenome((RSet))
predictSex <- getSex(GRSet, cutoff = -2)
predictSex$title <- rownames(predictSex)
predictSex <- as.data.frame(predictSex)

pheno <- left_join(pheno, predictSex, by = "title")
pheno$sex[pheno$sex == "male"] <- "M"
pheno$sex[pheno$sex == "female"] <- "F"
pheno <- mutate(pheno, sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally() #1 sample does not match
pheno$Sample_Name <- pheno$title

filter <- champ.filter(beta = getBeta(RSet),
                       pd = pheno,
                       Meth = meth,
                       UnMeth = unmeth,
                       detP = detP,
                       arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)


champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

#remove sample that sex does not match 

pheno <- pheno %>% filter(sexmatch == "Yes")
beta <- beta[,pheno$Sample_Name]

which(is.na(beta))

#normalisation 

myNorm <- champ.norm(beta = beta)

champ.QC(beta = myNorm,
         pheno = pheno$disease.status,
         dendrogram = FALSE)

champ.SVD(beta = myNorm,
          pd = pheno %>% select(age,
                                disease.status,
                                `levodopastatus:ch1`,
                                `numberofyearsinschool:ch1`,
                                `packyearstotal:ch1`,
                                ethnicity,
                                sex),
          resultsDir = "./ChAMP_SVDimages")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE72774 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

write.table(pheno,
            "GSE72774 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myNorm,
            "GSE72774 beta after normalisation filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")      


#### GSE72776 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE72776"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE72776",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
phenoraw <-pData(gset)
glimpse(phenoraw)

pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             `age:ch1`,
                             `b cell:ch1`,
                             `bioage4hastatic:ch1`,
                             `cd4.naive:ch1`,
                             `cd8.naive:ch1`,
                             `cd4t:ch1`,
                             `cd8.naive:ch1`,
                             `cd8pcd28ncd45ran:ch1`,
                             `cd8t:ch1`,
                             `disease status:ch1`,
                             `dnamage:ch1`,
                             `ethnicity:ch1`,
                             `granulocyte:ch1`,
                             `monocyte:ch1`,
                             `natural killer cell:ch1`,
                             `pdstudylevodopastatus:ch1`,
                             `plasmablast:ch1`,
                             `Sex:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "b cell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno) == "granulocyte:ch1"] <- "Gran"
names(pheno)[names(pheno) == "monocyte:ch1"] <- "Mono"
names(pheno)[names(pheno) == "plasmablast:ch1"] <- "Plasmablast"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "natural killer cell:ch1"] <- "NK"
names(pheno)[names(pheno) == "disease status:ch1"] <- "disease.status"
pheno$title <- sub("PEGhispanic: genomic DNA from whole blood of subject ","", pheno$title)

getGEOSuppFiles("GSE72776",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE72776/GSE72776"
setwd(path)

gunzip("GSE72776_datSignal.csv.gz")
gunzip("GSE72776_datBetaNormalized.csv.gz")

beta_processed <- read_csv("GSE72776_datBetaNormalized.csv")
signals <- read_csv("GSE72776_datSignal.csv")

#preprocess using signal intensities

rownames(signals) <- signals$TargetID
signals <- as.data.frame(signals)
meth <- signals %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Unmethylated","",colnames(meth))
meth <- meth %>% select(ends_with("Methylated"))
colnames(meth) <- sub("Methylated","", colnames(meth))
unmeth <- signals %>% select(contains("Unmethylated"))
colnames(unmeth) <- sub("Unmethylated","", colnames(unmeth))
detP <- signals %>% select(contains("DectionPvalue"))
colnames(detP) <- sub("DectionPvalue","", colnames(detP))
meth <- as.matrix(meth)
unmeth <- as.matrix(unmeth)
detP <- as.matrix(detP)

annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRSet <- mapToGenome((RSet))
predictSex <- getSex(GRSet, cutoff = -2)
predictSex$title <- rownames(predictSex)
predictSex <- as.data.frame(predictSex)

pheno <- left_join(pheno, predictSex, by = "title")
pheno$sex[pheno$sex == "male"] <- "M"
pheno$sex[pheno$sex == "female"] <- "F"
pheno <- mutate(pheno, sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally() #all samples match
pheno$Sample_Name <- pheno$title

filter <- champ.filter(beta = getBeta(RSet),
                       pd = pheno,
                       Meth = meth,
                       UnMeth = unmeth,
                       detP = detP,
                       arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)


champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

#normalisation 

myNorm <- champ.norm(beta = beta)

champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta = myNorm,
          pd = pheno %>% select(age,
                                disease.status,
                                ethnicity,
                                sex),
          resultsDir = "./ChAMP_SVDimages")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE72776 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

write.table(pheno,
            "GSE72776 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myNorm,
            "GSE72776 beta after normalisation filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")      


#### GSE166611 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE166611"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE166611",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)
pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          supplementary_file,
                          `ac:ch1`,
                          `age:ch1`,
                          `bmi:ch1`,
                          `cell type:ch1`,
                          `fm:ch1`,
                          `fm%:ch1`,
                          `gender:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "bmi:ch1"] <- "bmi"
names(pheno)[names(pheno) == "cell type:ch1"] <- "cell type"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"

pheno <- separate(pheno, col = supplementary_file, into = c("Supp info", "Sentrix_ID", "Sentrix_Position", "Channel colour"), sep = "\\_")

glimpse(pheno)
pheno$Sample_Name <- pheno$geo_accession
pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)
pheno$bmi <- as.numeric(pheno$bmi)

#download raw IDATS 

getGEOSuppFiles("GSE166611",
                fetch_files = TRUE)
untar("GSE166611/GSE166611_RAW.tar")

write.csv(pheno,
          "GSE166611 Phenotypes.csv")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))


filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$age <- as.numeric(pheno$age)

pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)
myNorm <- champ.norm(beta=beta,
                     arraytype="450K")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            bmi),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            bmi),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE166611 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE166611 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE166611 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE164056 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE164056"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE164056",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          supplementary_file,
                          `age:ch1`,
                          `sad:ch1`,
                          `ela:ch1`,
                          `ctq_total:ch1`,
                          `Sex:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "ctq_total:ch1"] <- "ctq_total"
names(pheno)[names(pheno) == "ela:ch1"] <- "early life adversity"
names(pheno)[names(pheno) == "sad:ch1"] <- "social anxiety disoder"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"

pheno$Sample_Name <- pheno$geo_accession
pheno <- separate(pheno, col = title, into = c("tissue", "Sentrix_ID", "Sentrix_Position"), sep = "\\_")

glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)


#download raw IDATS 

getGEOSuppFiles("GSE164056",
                fetch_files = TRUE)
untar("GSE164056/GSE164056_RAW.tar")

write.csv(pheno,
          "GSE164056 Phenotypes.csv")

path = "/pvol/Preprocessing/Blood/GEO/GSE164056"
setwd(path)
targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$age <- as.numeric(pheno$age)

predictedSex$Sample_Name <- rownames(predictedSex)
predictedSex <- as.data.frame(predictedSex)
pheno$Sample_Name <- rownames(pheno)
pheno <- pheno %>% left_join(predictedSex, by = "Sample_Name")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)
myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            social.anxiety.disoder,
                            early.life.adversity),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)
pheno$Slide <- as.factor(pheno$Slide)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            social.anxiety.disoder,
                            early.life.adversity),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

myCombat = logit2(myCombat)
pheno$Array <- as.factor(pheno$Array)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Array,
                mod=NULL)

#Run SVD again
myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            social.anxiety.disoder,
                            early.life.adversity),
          resultsDir="./CHAMP_SVDimages/batch_corrected/position corrected")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE164056 age distribution.tiff")

champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)


pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE164056 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE164056 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE85311 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE85311"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE85311",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

phenoraw = pData(gset)
glimpse(pheno)

pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             source_name_ch1,
                             supplementary_file,
                             description, 
                             description.1,
                             `age:ch1`,
                             `Sex:ch1`,
                             `training status:ch1`)
names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "training status:ch1"] <- "training status"

glimpse(pheno)

pheno$description.1 <- sub("YS","Ysed",pheno$description.1)
pheno$description.1 <- sub("OS","Osed",pheno$description.1)
pheno$description.1 <- sub("OT","Oex",pheno$description.1)


getGEOSuppFiles("GSE85311",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE85311/GSE85311"
setwd(path)

untar("GSE85311_RAW.tar")
gunzip("GSE85311_Unmethylated_and_methylated_signal_intensities.csv.gz")

rawsubset <- read.csv("GSE85311_Unmethylated_and_methylated_signal_intensities.csv", nrow = 5)
signals <- read.csv("GSE85311_Unmethylated_and_methylated_signal_intensities.csv")
rownames(signals) <- signals$X
meth <- signals %>% select(contains(".Methylated.Signal"))
meth <- as.matrix(meth)
colnames(meth) <- sub(".Methylated.Signal","",colnames(meth))
unmeth <- signals %>% select(contains(".Unmethylated.Signal"))
unmeth <- as.matrix(unmeth)
colnames(unmeth) <- sub(".Unmethylated.Signal","",colnames(unmeth))
detP <- signals %>% select(contains(".Detection.Pval"))
detP <- as.matrix(detP)
colnames(detP) <- sub(".Detection.Pval","",colnames(detP))
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictSex)
beta = getBeta(RSet)
beta <- as.matrix(beta)
pheno$Sample_Name <- pheno$description.1


filter <- champ.filter(beta = beta,
                       pd = pheno,
                       Meth = as.matrix(meth), 
                       UnMeth = as.matrix(unmeth), 
                       detP = as.matrix(detP),
                       arraytype = "450K")

champ.QC(beta = filter$beta,
         pheno = filter$pd$sex,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)

#check sexes

predictSex <- as.data.frame(predictSex)
predictSex$Sample_Name <- rownames(predictSex)
pheno <- left_join(pheno, predictSex, by = "Sample_Name")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

#normalisation of probes 

myNorm <- champ.norm(beta = beta, arraytype = "450K")

champ.QC(myNorm, 
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(as.matrix(myNorm), 
          pd = pheno %>% select(age,
                                sex,
                                training.status),
          resultsDir = "./ChaMP_SVDimages")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE85311 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

path = "/pvol/Preprocessing/Blood/GEO/GSE85311"
setwd(path)

write.table(pheno,
            "GSE85311 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myNorm,
            "GSE85311 beta after normalisation.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = TRUE)


#### GSE151278 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE151278"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE151278",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)
pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          supplementary_file,
                          `age at initiation:ch1`,
                          `age:ch1`,
                          `biological drug:ch1`,
                          `disease state:ch1`,
                          `gender:ch1`,
                          `response:ch1`,
                          `tissue:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "biological drug:ch1"] <- "drug"
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease state"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"

pheno <- separate(pheno, col = supplementary_file, into = c("Supp info", "Sentrix_ID", "Sentrix_Position", "Channel colour"), sep = "\\_")

glimpse(pheno)
pheno$age <- as.numeric(pheno$age)
pheno$`disease state` <- as.factor(pheno$`disease state`)
pheno$drug <- as.factor(pheno$drug)
pheno$Sample_Name <- paste0(pheno$geo_accession,"_",pheno$Sentrix_ID,"_",pheno$Sentrix_Position)

#download raw IDATS 

getGEOSuppFiles("GSE151278",
                fetch_files = TRUE,
                filter_regex = "tar$")

untar("GSE151278/GSE151278_RAW.tar")
getwd()
write.csv(pheno,
          "GSE151278 Phenotypes.csv")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))


filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"
predictedSex <- as.data.frame(predictedSex)
predictedSex$Sample_Name <- rownames(predictedSex)
pheno <- left_join(pheno, predictedSex, by = "Sample_Name")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="450K")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            drug,
                            `disease state`),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Sentrix_ID),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            drug,
                            `disease state`),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

M <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Sentrix_Position),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            drug,
                            `disease state`),
          resultsDir="./CHAMP_SVDimages/batch_Position_corrected/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE151278 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE151278 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE151278 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE96879 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE96879"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE96879",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)
pheno <- pheno %>% select(title, 
                          geo_accession,
                          supplementary_file,
                          `age:ch1`,
                          `ethnicity:ch1`,
                          `disease state:ch1`,
                          `sledai score:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "ethnicity:ch1"] <- "ethnicity"
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease state"
names(pheno)[names(pheno) == "sledai score:ch1"] <- "sledai score"

pheno <- separate(pheno, col = supplementary_file, into = c("Supp info", "Sentrix_ID", "Sentrix_Position", "Channel colour"), sep = "\\_")

glimpse(pheno)
pheno$age <- as.numeric(pheno$age)
pheno$`disease state` <- as.factor(pheno$`disease state`)
pheno$ethnicity<- as.factor(pheno$ethnicity)
pheno$`sledai score` <- as.double(pheno$`sledai score`)
pheno$Sample_Name <- paste0(pheno$geo_accession,"_",pheno$Sentrix_ID,"_",pheno$Sentrix_Position)

#download raw IDATS 

getGEOSuppFiles("GSE96879",
                fetch_files = TRUE,
                filter_regex = "tar$")

untar("GSE96879/GSE96879_RAW.tar")
getwd()
write.csv(pheno,
          "GSE96879 Phenotypes.csv")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

champ.QC(beta = beta,
         pheno = pheno$disease.state,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="450K")

champ.QC(beta = myNorm,
         pheno = pheno$disease.state,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            ethnicity,
                            disease.state,
                            sledai.score),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            ethnicity,
                            disease.state,
                            sledai.score),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")


pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE96879 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE96879 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE96879 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE134429 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE134429"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE134429",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

phenoraw = pData(gset)
glimpse(phenoraw)

pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             source_name_ch1,
                             characteristics_ch1,
                             description, 
                             description.1,
                             `age:ch1`,
                             `Sex:ch1`,
                             `batch:ch1`,
                             `donor:ch1`,
                             `das28:ch1`,
                             `patient cohort:ch1`)
names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "batch:ch1"] <- "batch"
names(pheno)[names(pheno) == "donor:ch1"] <- "donor"
names(pheno)[names(pheno) == "patient cohort:ch1"] <- "patient cohort"
names(pheno)[names(pheno) == "das28:ch1"] <- "das28"


getGEOSuppFiles("GSE134429",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE134429/GSE134429"
setwd(path)

gunzip("GSE134429_signal_intensities.txt.gz")

rawsubset <- read.delim("GSE134429_signal_intensities.txt", nrow = 5)
signals <- read.delim("GSE134429_signal_intensities.txt")
rownames(signals) <- signals$ID_REF
meth <- signals %>% select(contains("_Methylated_signal"))
meth <- as.matrix(meth)
colnames(meth) <- sub("_Methylated_signal","",colnames(meth))
unmeth <- signals %>% select(contains("_Unmethylated_signal"))
colnames(unmeth) <- sub("_Unmethylated_signal","",colnames(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% select(contains("_Detection_Pval"))
colnames(detP) <- sub("_Detection_Pval","",colnames(detP))
detP <- as.matrix(detP)
annotation <- c("IlluminaHumanMethylationEPIC","ilm10b4.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
pData(GRset)
predictSex <- getSex(GRset, cutoff = -2)
beta = getBeta(RSet)
beta <- as.matrix(beta)

pheno$Sample_Name <- sub("monocytes.","",pheno$title)
colnames(beta) == pheno$Sample_Name
colnames(beta)
pheno$Sample_Name
intersect(pheno$Sample_Name, colnames(detP))

filter <- champ.filter(beta = beta,
                       pd = pData(GRset),
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP,
                       arraytype = "EPIC")

champ.QC(beta = filter$beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


#check sexes

predictSex <- as.data.frame(predictSex)
predictSex$Sample_Name <- rownames(predictSex)
pheno <- left_join(pheno, predictSex, by = "Sample_Name")
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno <- pheno %>% filter(sexmatch == "Yes")
beta <- beta[,pheno$Sample_Name]
#normalisation of probes 

myNorm <- champ.norm(beta = beta, arraytype = "450K")

champ.QC(myNorm, 
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(myNorm, 
          pd = pheno %>% select(age,
                                sex,
                                batch,
                                das28,
                                `patient cohort`,
                                donor),
          resultsDir = "./ChaMP_SVDimages")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE134429 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

path = "/pvol/Preprocessing/Blood/GEO/GSE134429"
setwd(path)

write.table(pheno,
            "GSE134429 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myNorm,
            "GSE134429 beta after normalisation.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = TRUE)


#### GSE120307 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE120307"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE120307",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL13534	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
beta <- exprs(gset)
plotDensities(beta, legend = F)
phenoraw = pData(gset)
glimpse(pheno)
pheno <- phenoraw %>% select(title, 
                             geo_accession,
                             source_name_ch1,
                             supplementary_file,
                             `age:ch1`,
                             `diagnosis:ch1`,
                             `gender:ch1`,
                             `pair number:ch1`,
                             `tissue:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "diagnosis:ch1"] <- "diagnosis"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "pair number:ch1"] <- "twin pair"
names(pheno)[names(pheno) == "tissue:ch1"] <- "tissue"

glimpse(pheno)
pheno$age <- as.numeric(pheno$age)
pheno$Sample_Name <- pheno$geo_accession

#check for sexes 

data(probe.features)
x.probes <- filter(probe.features, CHR == "X")
x.probes <- rownames(x.probes)
predictSex <- wateRmelon::predictSex(beta, x.probes = x.probes, pc = 2, plot = TRUE)
predictSex <- as.data.frame(predictSex)
predictSex$Sample_Name <- rownames(predictSex)
pheno <- left_join(pheno, predictSex, by = "Sample_Name")
pheno$sex[pheno$sex == "male"] <- "Male"
pheno$sex[pheno$sex == "female"] <- "Female"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

filter = champ.filter(beta = beta,
                      pd =  pheno, 
                      arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


pheno = as.data.frame(filter$pd)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

impute <- champ.impute(beta = beta,
                       pd = pheno)

champ.QC(beta = impute$beta,
         pheno = impute$pd$sex,
         dendrogram = FALSE)

beta = impute$beta
pheno = impute$pd

champ.SVD(beta=beta,
          pd=pheno%>%select(age,
                            diagnosis,
                            sex,
                            `twin pair`),
          resultsDir="./CHAMP_SVDimages/")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE120307 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE120307 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(beta,
            "GSE120307 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE20236 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE20236"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE20236", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8490", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#extract phenotypes 

pheno <- pData(gset)
exprs <- exprs(gset)
plotDensities(exprs, legend = F)

#get raw supplementary files 

getGEOSuppFiles("GSE20236",
                fetch_files = TRUE)
gunzip("GSE20236/GSE20236_non-normalized.txt.gz")

path = "/pvol/Preprocessing/Blood/GEO/GSE20236/GSE20236"
setwd(path)
raw <- read_delim("GSE20236_non-normalized.txt")

#do not use raw files

#format pheno file 



names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
pheno <- pheno %>% select(title,
                          geo_accession,
                          age,
                          sex,
                          `tissue:ch1`)

pheno$Sample_Name <- pheno$geo_accession

filter <- champ.filter(beta = exprs,
                       pd = pheno, 
                       arraytype = "27K")

champ.QC(beta = filter$beta,
         pheno = filter$pd$sex,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}
pheno$age <- as.numeric(pheno$age)
pheno$age <- sub(" years","", pheno$age)
pheno %>% summarise(mean = mean(age),
                    sd = sd(age))

write.table(pheno, 
            "GSE20236 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")
write.table(beta,
            "GSE20236 beta after filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#perform champ.svd

B <- data.table::fread("GSE20236/GSE20236 beta after filtering.txt")
B <- as.data.frame(B)
rownames(B) <- B$V1
B <- B[,-1]
pheno <- read.delim("GSE20236/GSE20236 Phenotypes.txt")

champ.SVD(beta = as.matrix(B),
          pd = pheno)


#### GSE19711 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE19711"
setwd(path)

library(tidyverse)
library(ChAMP)
library(GEOquery)
library(limma)
library(minfi)

beta <- read.delim("GSE19711_beta.txt", sep = "", header = T)
pheno <- read.delim("GSE19711 Phenotypes.txt", sep = " ", header = T)
pheno_normal <- pheno %>%  subset(histology.ch1 == "")

keep <- pheno_normal$ID

beta_subset <- beta[,(names(beta) %in% keep)]

drop <- c("GSM492149","GSM492308")

beta_subset2 <- beta_subset[,!(names(beta_subset) %in% drop)]

plotDensities(beta_subset2, legend = F)

hist(beta_subset[,143])

colnames(beta_subset[11])
colnames(beta_subset[24])
colnames(beta_subset[25])
colnames(beta_subset[54])
colnames(beta_subset[61]) 
colnames(beta_subset[70]) 
colnames(beta_subset[143]) 

hist(beta_subset[,86])

glimpse(pheno)

write.table(beta_subset2, "GSE19711 beta outliers removed.txt")

pheno_normal <- pheno_normal %>% filter(ID != "GSM492149")
pheno_normal <- pheno_normal %>% filter(ID != "GSM492308")

glimpse(pheno_normal)

names(pheno_normal)[names(pheno_normal)=="age..at.recruitment."] <- "age"
names(pheno_normal)[names(pheno_normal)=="batch.ch1"] <- "batch"

glimpse(pheno_normal)

write.table(pheno_normal, "GSE19711 Phenotypes normal samples.txt")

#download raw files and preprocess manually 

getGEOSuppFiles("GSE19711",
                fetch_files = TRUE)

path = "/pvol/Preprocessing/Blood/GEO/GSE19711/GSE19711"
setwd(path)

untar("GSE19711_RAW.tar")
files <- list.files()[grep(".txt.gz$",list.files())]
for (f in files){
  gunzip(f)
}

raw <- list.files()[grep(".txt",list.files())]
signals <- data.frame(matrix("", nrow = 27578))  

for (r in raw){
  dat <- read.delim(r, header = TRUE, sep = "\t")
  r <- sub("\\.txt", "", r) #Obtain the name of the dataset without the ".txt" at the end
  print(r)
  colnames(dat) <- sub("",paste0(r,"."),colnames(dat))
  signals <- cbind(signals, dat)
}

rownames(signals) <- signals$GSM491937.IlmnID

beta <- signals %>% select(contains("Beta"))
rownames(beta)
names(beta) <- sub(".Beta","",names(beta))
meth <- signals %>% select(contains(".Methylated.channel.intensity"))
names(meth) <- sub(".Methylated.channel.intensity","",names(meth))
rownames(meth)
unmeth <- signals %>% select(contains(".Unmethylated.channel.intensity"))
names(unmeth) <- sub(".Unmethylated.channel.intensity","",names(unmeth))
rownames(unmeth)
detP <- signals %>% select(contains(".Detection.P.value"))
names(detP) <- sub(".Detection.P.value","",names(detP))
rownames(detP)

pheno_normal <- pheno %>% subset(histology.ch1 == "")
keep <- rownames(pheno_normal)
beta <- beta[,keep]
meth <- meth[,keep]
unmeth <- unmeth[,keep]

methylset <- MethylSet

filter <- champ.filter(beta = as.matrix(beta),
                       Meth = as.matrix(meth),
                       UnMeth = as.matrix(unmeth),
                       pd = pheno_normal,
                       detP = as.matrix(detP),
                       arraytype = "27K")
champ.QC(beta = filter$beta,
         pheno = filter$pd$ID,
         dendrogram = FALSE)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno <- filter$pd

impute <- champ.impute(beta,
                       pd = pheno)
beta <- impute$beta
beta <- as.matrix(beta)
pheno <- impute$pd

pheno <- separate(pheno, col = beadchip_well.ch1, into = c("Slide","Array"), sep = "\\_")

pheno <- pheno %>% select(ID, 
                          ageatdiagnosis.ch1,
                          age..at.recruitment.,
                          agegroupatsampledraw.ch1,
                          batch.ch1,
                          Slide,
                          Array,
                          bs.conversion.c1.ch1,
                          bs.conversion.c2.ch1,
                          ca125.ch1,
                          grade.ch1,
                          post.treatment.sample.ch1,
                          pre.treatment.sample.ch1,
                          sample.type.ch1,
                          stage.ch1)
names(pheno)[names(pheno) == "age..at.recruitment."] <- "age"
names(pheno)[names(pheno) == "ageatdiagnosis.ch1"] <- "age at diagnosis"
names(pheno)[names(pheno) == "agegroupatsampledraw.ch1"] <- "age group"
names(pheno)[names(pheno) == "batch.ch1"] <- "batch"

champ.SVD(beta = beta,
          pd = pheno %>% select(Slide,
                                Array,
                                age, 
                                batch,
                                `ca125.ch1`))
library(sva)
M <- logit2(beta)
pheno$Slide <- as.factor(pheno$Slide)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age, 
                            batch,
                            `ca125.ch1`),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

myCombat = logit2(myCombat)
pheno$Array <- as.factor(pheno$Array)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Array,
                mod=NULL)

#Run SVD again
myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age, 
                            batch,
                            `ca125.ch1`),
          resultsDir="./CHAMP_SVDimages/batch_corrected/position corrected")

write.table(pheno,
            "GSE19711 Phenotypes after preprocessing.txt",
            col.names = TRUE,
            sep = "\t")
write.table(myCombat,
            "GSE19711 beta after filtering and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### GSE157131 ####
#### GSE117859 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE117859"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE117859",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL16304",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          `age:ch1`,
                          `artadherence:ch1`,
                          `bcell:ch1`,
                          `cd4t:ch1`,
                          `cd8t:ch1`,
                          `gran:ch1`,
                          `hiv:ch1`,
                          `mono:ch1`,
                          `nk:ch1`,
                          `race:ch1`,
                          `Sex:ch1`,
                          `smoking:ch1`,
                          `wbc_new:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age" #converts samples without accurate ages into NAs. Remove from downstream analysis
names(pheno)[names(pheno) == "bcell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "gran:ch1"] <- "Gran"
names(pheno)[names(pheno) == "hiv:ch1"] <- "HIV"
names(pheno)[names(pheno) == "mono:ch1"] <- "Mono"
names(pheno)[names(pheno) == "nk:ch1"] <- "NK"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "smoking:ch1"] <- "smoking"
names(pheno)[names(pheno) == "race:ch1"] <- "race"

glimpse(pheno)

pheno$Sample_Name <- pheno$geo_accession
pheno <- separate(pheno, col = source_name_ch1, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")

glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)
pheno <- pheno %>% modify_at(c("CD4T","CD8T","Bcell","Mono","Gran","NK"), as.numeric)


#download raw files 

getGEOSuppFiles("GSE117859",
                fetch_files = TRUE)
gunzip("GSE117859/GSE117859_MethylatedSignal.txt.gz")

signals <- read.delim("GSE117859/GSE117859_MethylatedSignal.txt")
rownames(signals) <- signals$ID_REF
meth <- signals %>% select(contains(".Methylated.signal"))
colnames(meth) <- sub(".Methylated.signal","",colnames(meth))
meth <- as.matrix(meth)
unmeth <- signals %>% select(contains(".Unmethylated.signal"))
colnames(unmeth) <- sub(".Unmethylated.signal","",colnames(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% select(contains(".Detection.Pval"))
colnames(detP) <- sub(".Detection.Pval","",colnames(detP))
detP <- as.matrix(detP)
library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
getSex(GRset, cutoff = -2)
beta = getBeta(RSet)
beta <- as.matrix(beta)

pheno <- separate(pheno, col = title, into = c("Sample_Name", "addition info","add info 2"), sep = "\\_")
pheno$Sample_Name <- sub(" ",".", pheno$Sample_Name)

filter <- champ.filter(beta = beta,
                       pd = pheno,
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP,
                       arraytype = "450K")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)


champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="450K")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)
myNorm <- as.matrix(myNorm)
champ.SVD(beta=myNorm,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)
pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/batch_corrected")

myCombat <- logit2(myCombat)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

#Run SVD again
myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/batch_corrected/position corrected")


pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE147740 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

write.table(pheno,
            "GSE117859 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE117859 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### GSE117860 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE117860"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE117860",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          `age:ch1`,
                          `artadherence:ch1`,
                          `bcell_850k:ch1`,
                          `cd4t_850k:ch1`,
                          `cd8t_850k:ch1`,
                          `gran_850k:ch1`,
                          `hiv:ch1`,
                          `mono_850k:ch1`,
                          `nk_850k:ch1`,
                          `race:ch1`,
                          `Sex:ch1`,
                          `smoking:ch1`,
                          `wbc_new:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age" #converts samples without accurate ages into NAs. Remove from downstream analysis
names(pheno)[names(pheno) == "bcell_850k:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t_850k:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t_850k:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "gran_850k:ch1"] <- "Gran"
names(pheno)[names(pheno) == "hiv:ch1"] <- "HIV"
names(pheno)[names(pheno) == "mono_850k:ch1"] <- "Mono"
names(pheno)[names(pheno) == "nk_850k:ch1"] <- "NK"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "smoking:ch1"] <- "smoking"
names(pheno)[names(pheno) == "race:ch1"] <- "race"

glimpse(pheno)

pheno <- separate(pheno, col = source_name_ch1, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")

glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)
pheno <- pheno %>% modify_at(c("CD4T","CD8T","Bcell","Mono","Gran","NK"), as.numeric)


#download raw files 

getGEOSuppFiles("GSE117860",
                fetch_files = TRUE)
gunzip("GSE117860/GSE117860_MethylatedSignal.txt.gz")

signals <- read.delim("GSE117860/GSE117860_MethylatedSignal.txt")
rownames(signals) <- signals$ID_REF
meth <- signals %>% select(contains(".Methylated.signal"))
colnames(meth) <- sub(".Methylated.signal","",colnames(meth))
meth <- as.matrix(meth)
unmeth <- signals %>% select(contains(".Unmethylated.signal"))
colnames(unmeth) <- sub(".Unmethylated.signal","",colnames(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% select(contains(".Detection.Pval"))
colnames(detP) <- sub(".Detection.Pval","",colnames(detP))
detP <- as.matrix(detP)
library(minfi)
annotation <- c("IlluminaHumanMethylationEPIC","ilm10b4.hg19")
names(annotation) <- c("array","annotation")
#check which cpgs have NAs in the detection p values
detP_NA <- apply(detP,1,function(x){any(is.na(x))})

methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
getSex(GRset, cutoff = -2)
beta = getBeta(RSet)
beta <- as.matrix(beta)

#remove the NAs from the detP values 
detP <- detP[!detP_NA,]
beta <- beta[!detP_NA,]
meth <- meth[!detP_NA,]
unmeth <- unmeth[!detP_NA,]
pheno <- pheno %>% filter(Sample_Name %in% colnames(beta))
pheno <- separate(pheno, col = title, into = c("Sample_Name", "addition info","add info 2"), sep = "\\_")
pheno$Sample_Name <- sub(" ",".", pheno$Sample_Name)

filter <- champ.filter(beta = beta,
                       pd = pheno,
                       Meth = meth, 
                       UnMeth = unmeth,
                       arraytype = "EPIC")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)


champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

#can't normalise the beta values?

impute <- champ.impute(beta = beta,
                       pd = pheno)

myNorm <- champ.norm(beta=impute$beta,arraytype="EPIC")

champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)
pheno <- impute$pd

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)
pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/batch_corrected")

myCombat <- logit2(myCombat)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

#Run SVD again
myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Sentrix_ID,
                            Sentrix_Position,
                            age,
                            sex,
                            smoking,
                            HIV),
          resultsDir="./CHAMP_SVDimages/batch_corrected/position corrected")


pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE147740 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))


champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

write.table(pheno,
            "GSE117860 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE117860 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE147740 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE147740"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE147740",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title, 
                          geo_accession,
                          source_name_ch1,
                          supplementary_file,
                          `age:ch1`,
                          `disease state:ch1`,
                          `person id:ch1`,
                          `Sex:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age" #converts samples without accurate ages into NAs. Remove from downstream analysis
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease.state"
names(pheno)[names(pheno) == "person id:ch1"] <- "person.id"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
glimpse(pheno)

pheno$Sample_Name <- pheno$geo_accession
pheno <- separate(pheno, col = supplementary_file, into = c("Supp_info", "Sentrix_ID", "Sentrix_Position","colour channel"), sep = "\\_")

glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno$sex <- as.factor(pheno$sex)

#can't identify replicates with phenotype information
phenoraw = pData(gset)
glimpse(phenoraw)


#download raw IDATS 

getGEOSuppFiles("GSE147740",
                fetch_files = TRUE)
untar("GSE147740/GSE147740_RAW.tar")

write.csv(pheno,
          "GSE147740 Phenotypes.csv")

path = "/pvol/Preprocessing/Blood/GEO/GSE147740"
setwd(path)
targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")
library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)

pheno$age <- as.numeric(pheno$age)
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno <- pheno %>% filter(sexmatch == "Yes")
toKeep <- rownames(pheno)
beta <- beta[,toKeep]

pheno <- pheno %>% filter(!is.na(age))
toKeep <- rownames(pheno)
beta <- beta[,toKeep]

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)
myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state),
          resultsDir="./CHAMP_SVDimages/")
library(sva)
M <- logit2(myNorm)
pheno$Slide <- as.factor(pheno$Slide)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Slide,
                mod=NULL)

myCombat=ilogit2(myCombat)

#Run SVD again
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state),
          resultsDir="./CHAMP_SVDimages/batch_corrected/")

myCombat = logit2(myCombat)
pheno$Array <- as.factor(pheno$Array)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Array,
                mod=NULL)

#Run SVD again
myCombat = ilogit2(myCombat)
champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state),
          resultsDir="./CHAMP_SVDimages/batch_corrected/position corrected")

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE147740 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
head(pheno)[,1:5]
write.table(pheno,
            "GSE147740 Phenotypes.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE147740 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#average the replicates 
unique <- unique(pheno$person.id)
replicates <- pheno %>% group_by(person.id) %>% summarise(n = n()) %>% filter(n>1) %>% arrange(-n)
replicate_ids <- replicates$person.id
pheno$geo_accession <- rownames(pheno)

gsms <- c()
pheno_toadd <- rbind()
avgs <- cbind()
for (d in replicate_ids){
  subgsms <- pheno %>% filter(person.id == d) %>% pull(geo_accession)
  gsms <- c(gsms, subgsms)
  pheno_toadd <- rbind(pheno_toadd,
                       pheno %>% filter(geo_accession==subgsms[1]))
  avgs <- cbind(avgs,rowMeans(myCombat[,subgsms]))
}

colnames(avgs) <- replicate_ids
myCombat_avg <- myCombat
pheno_avg <- pheno

myCombat_avg <- cbind(myCombat_avg,
                      avgs)

pheno_avg <- bind_rows(pheno_avg,
                       pheno_toadd)
gsms_to_change <- tibble(geo_accession = colnames(myCombat_avg)[grepl("GSM",colnames(myCombat_avg))])
gsms_to_change <- left_join(gsms_to_change, pheno) %>% pull(person.id)
colnames(myCombat_avg)[grepl("GSM",colnames(myCombat_avg))] <- gsms_to_change

write.table(myCombat_avg,
            "GSE147740 normalised and batch corrected replicates avg beta.txt",
            col.names = TRUE,
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")

pheno_avg$person.id == colnames(myCombat_avg)

write.table(pheno_avg,
            "GSE147740 Phenotypes with average replicates.txt",
            quote = FALSE,
            col.names = TRUE,
            sep = "\t")

pheno %>% group_by(sex) %>% tally()


#### GSE152026 ####
path = "/pvol/Preprocessing/Blood/GEO/GSE152026"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)

gset <- getGEO("GSE152026",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

getGEOSuppFiles("GSE152026",
                fetch_files = TRUE,
                filter_regex = "GSE152026_EUGEI_raw_Signal.csv.gz$")

gunzip("GSE152026/GSE152026_EUGEI_raw_Signal.csv.gz")


pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          description,
                          `age:ch1`,
                          `phenotype:ch1`,
                          `Sex:ch1`)
pheno$Sample_Name <- pheno$title
pheno <- separate(pheno, col = Sample_Name, into = c("Sample_Name", "Samp_description"), sep = "\\ ")
pheno <- separate(pheno, col= title, into = c("Sentrix_ID","Sentrix_Position", sep = "\\ "))

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "phenotype:ch1"] <- "disease status"

getwd()

write.table(pheno,
            "GSE152025 Phenotypes.txt",
            sep = "\t")

signals <- data.table::fread("GSE152026/GSE152026_EUGEI_raw_Signal.csv", sep = ",")
rownames(signals) <- signals$V1
meth <- signals %>% select(contains("_Methylated_Signal"))
names(meth) <- sub("_Methylated_Signal","",names(meth))
meth <- as.matrix(meth)
rownames(meth) <- signals$V1
unmeth <- signals %>% select(contains("_Unmethylated_Signal"))
names(unmeth) <- sub("_Unmethylated_Signal","",names(unmeth))
unmeth <- as.matrix(unmeth)
rownames(unmeth) <- signals$V1
detP <- signals %>% select(contains("_Detection_Pval"))
names(detP) <- sub("_Detection_Pval","",names(detP))
detP <- as.matrix(detP)
rownames(detP) <- signals$V1

#check for sex
library(minfi)
annotation <- c("IlluminaHumanMethylationEPIC","ilm10b4.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset,
                       cutoff = -2)
betaraw = getBeta(RSet)
betaraw <- as.matrix(betaraw)

pheno <- read.delim("GSE152025 Phenotypes.txt")
pheno$Sample_Name == colnames(betaraw)
Sample_Name <- colnames(betaraw)
Sample_Name <- as.data.frame(Sample_Name)
pheno <- left_join(Sample_Name, pheno, by = "Sample_Name")

pheno$Sample_Name == colnames(betaraw)


filter <- champ.filter(beta = betaraw,
                       pd = as.data.frame(pheno),
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP,
                       arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

champ.QC(beta, 
         pheno = filter$pd$sex,
         dendrogram = FALSE)

predictedSex$Sample_Name <- rownames(predictedSex)
predictedSex <- as.data.frame(predictedSex)
pheno <- left_join(pheno, predictedSex, by = "Sample_Name")
pheno$sex[pheno$sex == "Female"] <- "F"
pheno$sex[pheno$sex == "Male"] <- "M"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno <- pheno %>% filter(sexmatch == "Yes")
keep <- pheno$Sample_Name
beta <- beta[,keep]


myNorm <- champ.norm(beta = beta, arraytype = "EPIC")

champ.QC(myNorm, 
         pheno = pheno$sex,
         dendrogram = FALSE)

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE152026 age distribution.tiff")

pheno %>% filter(!is.na(age)) %>% summarise(mean = mean(age),
                                            sd = sd(age),
                                            range = range(age))

pheno %>% group_by(sex) %>% tally()

champ.SVD(myNorm, 
          pd = pheno %>% select(age,
                                sex,
                                disease.status,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages")

pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)
pheno$age <- as.numeric(pheno$age)
pheno$disease.status <- as.factor(pheno$disease.status)
pheno$sex <- as.factor(pheno$sex)

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

myCombat=ilogit2(myCombat)


write.table(myCombat,
            file="GSE152026 beta normalised and batch corrected.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(pheno,
            file="GSE152026 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')


#### GSE132203 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE132203"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)
library(minfi)

gset <- getGEO("GSE132203",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          `age acceleration:ch1`,
                          `age:ch1`,
                          `bcell:ch1`,
                          `cd4t:ch1`,
                          `cd8t:ch1`,
                          `childabphyssexemot_ctq_01modandsev:ch1`,
                          `gender:ch1`,
                          `mergedcapsandpsswinthin30days:ch1`,
                          `mono:ch1`,
                          `neu:ch1`,
                          `nk:ch1`,
                          `pc1:ch1`,
                          `pc2:ch1`,
                          `pc3:ch1`,
                          `race:ch1`,
                          `tei_total_types_experienced_somewitness:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "bcell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "mono:ch1"] <- "Mono"
names(pheno)[names(pheno) == "neu:ch1"] <- "Neutro"
names(pheno)[names(pheno) == "nk:ch1"] <- "NK"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "race:ch1"] <- "ethnicity"

pheno$Sentrix_ID <- pheno$title
pheno <- separate(pheno, col = Sentrix_ID, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")
pheno$Sample_Name <- pheno$title


write.csv(pheno,
          "GSE132203 Phenotypes.csv")

getGEOSuppFiles("GSE132203",
                fetch_files = TRUE,
                filter_regex = ".tar$")
untar("GSE132203/GSE132203_RAW.tar")

getwd()


targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
glimpse(pheno)
pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE132203 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

glimpse(pheno)
champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            childabphyssexemot_ctq_01modandsev.ch1,
                            mergedcapsandpsswinthin30days.ch1,
                            tei_total_types_experienced_somewitness.ch1,
                            ethnicity),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            childabphyssexemot_ctq_01modandsev.ch1,
                            mergedcapsandpsswinthin30days.ch1,
                            tei_total_types_experienced_somewitness.ch1,
                            ethnicity),
          resultsDir="./CHAMP_SVDimages/")

myCombat=ilogit2(myCombat)

write.table(pheno,
            "GSE132203 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE132203 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE100264 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE100264"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE100264",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL16304",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

getGEOSuppFiles("GSE100264",
                fetch_files = TRUE,
                filter_regex = "MethylatedSignal.txt.gz$")

gunzip("GSE100264/GSE100264_MethylatedSignal.txt.gz")


pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          description,
                          `age:ch1`,
                          `artadherence:ch1`,
                          `bcell:ch1`,
                          `cd4t:ch1`,
                          `cd8t:ch1`,
                          `gran:ch1`,
                          `hcv_dx:ch1`,
                          `hiv:ch1`,
                          `idu:ch1`,
                          `mono:ch1`,
                          `nk:ch1`,
                          `race:ch1`,
                          `smoking:ch1`,
                          `Sex:ch1`)

pheno$Sample_Name <- pheno$title
pheno <- separate(pheno, col= title, into = c("Sentrix_ID","Sentrix_Position"), sep = "\\_")

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "bcell:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "gran:ch1"] <- "Gran"
names(pheno)[names(pheno) == "mono:ch1"] <- "Mono"
names(pheno)[names(pheno) == "nk:ch1"] <- "NK"
names(pheno)[names(pheno) == "smoking:ch1"] <- "Smoking"
names(pheno)[names(pheno) == "race:ch1"] <- "race"
names(pheno)[names(pheno) == "hiv:ch1"] <- "HIV"
names(pheno)[names(pheno) == "hcv_dx:ch1"] <- "hcv_dx"
names(pheno)[names(pheno) == "artadherence:ch1"] <- "artadherence"
names(pheno)[names(pheno) == "idu:ch1"] <- "idu"

getwd()

write.table(pheno,
            "GSE100264 Phenotypes.txt",
            sep = "\t")

pheno$Sample_Name <- pheno$description

signals <- data.table::fread("GSE100264/GSE100264_MethylatedSignal.txt")
signals <- as.data.frame(signals)
rownames(signals) <- signals$ID_REF
meth <- signals %>% select(contains(" Methylated signal"))
names(meth) <- sub(" Methylated signal","",names(meth))
meth <- as.matrix(meth)
unmeth <- signals %>% select(contains(" Unmethylated signal"))
names(unmeth) <- sub(" Unmethylated signal","",names(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% select(contains(" Detection Pval"))
names(detP) <- sub(" Detection Pval","",names(detP))
detP <- as.matrix(detP)

#check for sex
library(minfi)
annotation <- c("IlluminaHumanMethylation450k","ilmn12.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset,
                       cutoff = -2)
beta = getBeta(RSet)
beta <- as.matrix(beta)

pheno <- read.delim("GSE100264 Phenotypes.txt")
pheno$Sample_Name == colnames(beta)
#Sample_Name <- colnames(betaraw)
#Sample_Name <- as.data.frame(Sample_Name)
#pheno <- left_join(Sample_Name, pheno, by = "Sample_Name")

#pheno$Sample_Name == colnames(betaraw)


filter <- champ.filter(beta = beta,
                       pd = as.data.frame(pheno),
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP,
                       arraytype = "450K")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

champ.QC(beta, 
         pheno = filter$pd$sex,
         dendrogram = FALSE)

predictedSex$Sample_Name <- rownames(predictedSex)
predictedSex <- as.data.frame(predictedSex)
pheno <- left_join(pheno, predictedSex, by = "Sample_Name")
#pheno$sex[pheno$sex == "Female"] <- "F"
pheno$sex[pheno$sex == "Male"] <- "M"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

#pheno <- pheno %>% filter(sexmatch == "Yes")
#keep <- pheno$Sample_Name
#beta <- beta[,keep]

myNorm <- champ.norm(beta = beta, arraytype = "450K")

champ.QC(myNorm, 
         pheno = pheno$sex,
         dendrogram = FALSE)

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE100264 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

champ.SVD(myNorm, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages")

pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(myCombat, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages/batchcorrected")

myCombat <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(myCombat, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages/batchcorrected_positioncorrected")

write.table(myCombat,
            file="GSE100264 beta normalised and batch corrected.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(pheno,
            file="GSE100264 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')



#### GSE107080 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE107080"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)


gset <- getGEO("GSE107080",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

getGEOSuppFiles("GSE107080",
                fetch_files = TRUE,
                filter_regex = "MethylatedSignal.txt.gz$")

gunzip("GSE107080/GSE107080_MethylatedSignal.txt.gz")


pheno= pData(gset)

glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          description,
                          description.2,
                          `age:ch1`,
                          `artadherence:ch1`,
                          `bcell_850k:ch1`,
                          `cd4t_850k:ch1`,
                          `cd8t_850k:ch1`,
                          `gran_850k:ch1`,
                          `hcv_dx:ch1`,
                          `hiv:ch1`,
                          `idu:ch1`,
                          `mono_850k:ch1`,
                          `nk_850k:ch1`,
                          `race:ch1`,
                          `smoking:ch1`,
                          `Sex:ch1`,
                          `wbc_new:ch1`)

pheno <- separate(pheno, col= description, into = c("Sentrix_ID","Sentrix_Position","Colour"), sep = "\\_")

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "Sex:ch1"] <- "sex"
names(pheno)[names(pheno) == "bcell_850k:ch1"] <- "Bcell"
names(pheno)[names(pheno) == "cd4t_850k:ch1"] <- "CD4T"
names(pheno)[names(pheno) == "cd8t_850k:ch1"] <- "CD8T"
names(pheno)[names(pheno) == "gran_850k:ch1"] <- "Gran"
names(pheno)[names(pheno) == "mono_850k:ch1"] <- "Mono"
names(pheno)[names(pheno) == "nk_850k:ch1"] <- "NK"
names(pheno)[names(pheno) == "smoking:ch1"] <- "Smoking"
names(pheno)[names(pheno) == "race:ch1"] <- "race"
names(pheno)[names(pheno) == "hiv:ch1"] <- "HIV"
names(pheno)[names(pheno) == "hcv_dx:ch1"] <- "hcv_dx"
names(pheno)[names(pheno) == "artadherence:ch1"] <- "artadherence"
names(pheno)[names(pheno) == "idu:ch1"] <- "idu"

getwd()


pheno$Sample_Name <- pheno$description.2
write.table(pheno,
            "GSE107080 Phenotypes.txt",
            sep = "\t")

signals <- data.table::fread("GSE107080/GSE107080_MethylatedSignal.txt")
signals <- as.data.frame(signals)
rownames(signals) <- signals$ID_REF
meth <- signals %>% select(contains(" Methylated signal"))
names(meth) <- sub(" Methylated signal","",names(meth))
meth <- as.matrix(meth)
unmeth <- signals %>% select(contains(" Unmethylated signal"))
names(unmeth) <- sub(" Unmethylated signal","",names(unmeth))
unmeth <- as.matrix(unmeth)
detP <- signals %>% select(contains(" Detection Pval"))
names(detP) <- sub(" Detection Pval","",names(detP))
detP <- as.matrix(detP)

#check for sex
library(minfi)
annotation <- c("IlluminaHumanMethylationEPIC","ilm10b4.hg19")
names(annotation) <- c("array","annotation")
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset,
                       cutoff = -2)
beta = getBeta(RSet)
beta <- as.matrix(beta)

pheno <- read.delim("GSE107080 Phenotypes.txt")
pheno$Sample_Name == colnames(beta)
#Sample_Name <- colnames(betaraw)
#Sample_Name <- as.data.frame(Sample_Name)
#pheno <- left_join(Sample_Name, pheno, by = "Sample_Name")

#pheno$Sample_Name == colnames(betaraw)

row.has.na <- apply(beta,1,function(x){any(is.na(x))})
beta_filterednas <- beta[!row.has.na,]

#upload 450K dataset to check probes

Beta_GSE100264 <- data.table::fread("/pvol/Preprocessing/Blood/GEO/GSE100264/GSE100264 beta normalised and batch corrected.txt")
probes_450K <- Beta_GSE100264$V1
probes_450K <- as.data.frame(probes_450K)
probes_beta <- rownames(beta_filterednas)
probes_beta <- as.data.frame(probes_beta)
probes_beta$probes <- probes_beta$probes_beta
probes_450K$probes <- probes_450K$probes_450K
matchingprobes <- inner_join(probes_450K, probes_beta)

meth <- meth[!row.has.na,]
unmeth <- unmeth[!row.has.na,]
detP <- detP[!row.has.na,]

filter <- champ.filter(beta = beta_filterednas,
                       pd = pheno,
                       Meth = meth, 
                       UnMeth = unmeth, 
                       detP = detP)

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


champ.QC(beta, 
         pheno = pheno$sex,
         dendrogram = FALSE)

#predictedSex$Sample_Name <- rownames(predictedSex)
#predictedSex <- as.data.frame(predictedSex)
#pheno <- left_join(pheno, predictedSex, by = "Sample_Name")
#pheno$sex[pheno$sex == "Female"] <- "F"
#pheno$sex[pheno$sex == "Male"] <- "M"
#pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
#pheno %>% group_by(sexmatch) %>% tally()

#pheno <- pheno %>% filter(sexmatch == "Yes")
#keep <- pheno$Sample_Name
#beta <- beta[,keep]

myNorm <- champ.norm(beta = beta)

champ.QC(myNorm, 
         pheno = pheno$sex,
         dendrogram = FALSE)

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE107080 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))

pheno %>% group_by(sex) %>% tally()

champ.SVD(myNorm, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages")

pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
pheno$Sentrix_Position <- as.factor(pheno$Sentrix_Position)

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(myCombat, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages/batchcorrected")

myCombat <- logit2(myCombat)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.SVD(myCombat, 
          pd = pheno %>% select(age,
                                artadherence,
                                hcv_dx,
                                HIV,
                                idu,
                                Smoking,
                                race,
                                Sentrix_ID,
                                Sentrix_Position),
          resultsDir = "./CHAMP_SVDimages/batchcorrected_positioncorrected")

write.table(myCombat,
            file="GSE107080 beta normalised and batch corrected.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')

write.table(pheno,
            file="GSE107080 Phenotypes.txt",
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep='\t')



#### GSE116339 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE116339"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)
library(minfi)

gset <- getGEO("GSE116339",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          `age:ch1`,
                          `gender:ch1`,
                          `ln(totalpbb):ch1`,
                          `pbb-101:ch1`,
                          `pbb-153:ch1`,
                          `pbb-180:ch1`,
                          `pbb-77:ch1`,
                          `sample_id:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "sample_id:ch1"] <- "Sample_Name"

pheno <- separate(pheno, col = title, into = c("Sentrix_ID", "Sentrix_Position"), sep = "\\_")

write.csv(pheno,
          "GSE116339 Phenotypes.csv")

getGEOSuppFiles("GSE116339",
                fetch_files = TRUE,
                filter_regex = ".tar$")
untar("GSE116339/GSE116339_RAW.tar")

getwd()

targets <- read.metharray.sheet(getwd(),"GSE116339 Phenotypes.csv")
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
pheno$Sample_Name <- rownames(pheno)
pheno$sex[pheno$sex == "Male"] <- "M"
pheno$sex[pheno$sex == "Female"] <- "F"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno <- pheno %>% filter(sexmatch == "Yes")
keep <- pheno$Sample_Name
beta <- as.data.frame(beta)
beta <- beta[,keep]

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE116339 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

champ.QC(beta = as.matrix(beta),
         pheno = pheno$pbb.101.ch1,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")
champ.QC(beta = myNorm,
         pheno = pheno$pbb.180.ch1,
         dendrogram = FALSE)

glimpse(pheno)

#download their preprocessed matrix and compare
getGEOSuppFiles("GSE116339",
                fetch_files = TRUE,
                filter_regex = ".csv.gz$")

gunzip("GSE116339/GSE116339_Processed_Matrix.csv.gz")
preprocessed_beta <- data.table::fread("GSE116339/GSE116339_Processed_Matrix.csv", sep = ",")
preprocessed_beta <- as.data.frame(preprocessed_beta)
rownames(preprocessed_beta) <- preprocessed_beta$V1
preprocessed_beta <- select(preprocessed_beta, -V1)

limma::plotDensities(preprocessed_beta,legend = F)

#rather use their preprocessed beta matrix 

colnames(preprocessed_beta) <- sub("X","",colnames(preprocessed_beta))

champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            ln.totalpbb..ch1,
                            pbb.101.ch1,
                            pbb.153.ch1,
                            pbb.180.ch1,
                            pbb.77.ch1),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            ln.totalpbb..ch1,
                            pbb.101.ch1,
                            pbb.153.ch1,
                            pbb.180.ch1,
                            pbb.77.ch1),
          resultsDir="./CHAMP_SVDimages/batchandpositioncorrected")

myCombat=ilogit2(myCombat)


champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

write.table(pheno,
            "GSE116339 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE116339 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE168739 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE168739"
setwd(path)

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(limma)
library(minfi)


gset <- getGEO("GSE168739",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145	",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno = pData(gset)
glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          supplementary_file,
                          `age:ch1`,
                          `gender:ch1`,
                          `disease state:ch1`,
                          `individual id:ch1`,
                          `tissue:ch1`)

names(pheno)[names(pheno) == "age:ch1"] <- "age"
names(pheno)[names(pheno) == "gender:ch1"] <- "sex"
names(pheno)[names(pheno) == "individual id:ch1"] <- "individual id"
names(pheno)[names(pheno) == "disease state:ch1"] <- "disease.state"

pheno <- separate(pheno, col = supplementary_file, into = c("Suppinfo","Sentrix_ID", "Sentrix_Position","Colour"), sep = "\\_")

write.csv(pheno,
          "GSE168739 Phenotypes.csv")

getGEOSuppFiles("GSE168739",
                fetch_files = TRUE,
                filter_regex = ".tar$")
untar("GSE168739/GSE168739_RAW.tar")

getwd()

targets <- read.metharray.sheet(getwd(),"GSE168739 Phenotypes.csv")
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)
beta = getBeta(RSet)
pheno = as.matrix(pData(GRset))

filter = champ.filter(beta = beta,
                      pd =  pheno,
                      detP = detP,
                      Meth = meth,
                      UnMeth = unmeth, 
                      arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)
pheno$Sample_Name <- rownames(pheno)
#pheno$sex[pheno$sex == "Male"] <- "M"
#pheno$sex[pheno$sex == "Female"] <- "F"
pheno <- pheno %>% mutate(sexmatch = ifelse(sex == predictedSex,"Yes","No"))
pheno %>% group_by(sexmatch) %>% tally()

pheno <- pheno %>% filter(sexmatch == "Yes")
keep <- pheno$Sample_Name
beta <- as.data.frame(beta)
beta <- beta[,keep]

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE18739 age distribution.tiff")

pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

champ.QC(beta = as.matrix(beta),
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")
champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

glimpse(pheno)

champ.SVD(beta=as.matrix(myNorm),
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex,
                            disease.state),
          resultsDir="./CHAMP_SVDimages/batchandpositioncorrected")

myCombat=ilogit2(myCombat)


champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

write.table(pheno,
            "GSE168739 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE168739 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")


#### GSE197674 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE197674"
setwd(path)

library(tidyverse)
library(ChAMP)
library(limma)
library(minfi)
library(data.table)
library(GEOquery)

#### preprocessing 

gset <- getGEO("GSE197674",
               destdir = getwd(),
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL21145",  #this is the specific platform e.g. 450K or EPIC
                                  attr(gset, "names")) else idx <- 1
gset <- gset[[1]]

pheno = pData(gset)
glimpse(pheno)

getGEOSuppFiles("GSE197674",
                fetch_files = TRUE,
                filter_regex = "processed.txt.gz$")

gunzip("GSE197674/GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_processed.txt.gz")

beta <- data.table::fread("GSE197674/GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_processed.txt",
                          sep = "\t")
beta <- as.data.frame(beta)
rownames(beta) <- beta$ID_REF
beta <- select(beta,!contains("Detection Pval"))
rownames(beta)
beta <- select(beta,!contains("ID_REF"))

glimpse(pheno)

pheno <- pheno %>% select(title,
                          geo_accession,
                          supplementary_file,
                          `abdominal_pelvic_rt:ch1`,
                          `age:ch1`,
                          `alkylating agent,classic:ch1`,
                          `anthracyclines:ch1`,
                          `brainrt:ch1`,
                          `chestrt:ch1`,
                          `corticosteroids:ch1`,
                          `epipodophyllotoxins:ch1`,
                          `gender:ch1`,
                          `platinum:ch1`,
                          `vincristine:ch1`)

pheno <- separate(pheno, col = title, into = c("title","Sample_Name"), sep = "\\[")
pheno$Sample_Name <- sub("]","",pheno$Sample_Name)
names(pheno)[names(pheno)=="age:ch1"] <- "age"
names(pheno)[names(pheno)=="gender:ch1"] <- "sex"

pheno$Sample_Name == colnames(beta)
Sample_Name <- colnames(beta)
Sample_Name <- as.data.frame(Sample_Name)

pheno <- Sample_Name %>% left_join(pheno, by = "Sample_Name")
pheno$Sample_Name == colnames(beta)

gc()
champ.QC(beta = as.matrix(beta),
         pheno = pheno$sex,
         dendrogram = FALSE)

#predict sex using watermelon
library(wateRmelon)
data(probe.features)
Xprobes <- probe.features %>% filter(CHR == "X")
Xprobes <- rownames(Xprobes)
betaX <- beta[Xprobes,]
predictSex <- predictSex(as.matrix(betaX), pc = 2, plot = TRUE)
champ.QC()

write.table(pheno,
            "GSE197674 Phenotypes.txt")

filter <- champ.filter(beta = as.matrix(beta),
                       pd = pheno, 
                       arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno = as.data.frame(filter$pd)

pheno$age <- as.numeric(pheno$age)
ages_range <- ggplot(data = pheno) +
  geom_histogram(mapping = aes(x = age), colour = "black", fill = "skyblue", bins = 40) + 
  theme_minimal()
ages_range

ggsave("GSE197674 age distribution.tiff")


pheno %>% summarise(mean = mean(age),
                    sd = sd(age),
                    range = range(age))
pheno %>% group_by(sex) %>% tally()

champ.QC(beta = as.matrix(beta),
         pheno = pheno$sex,
         dendrogram = FALSE)

write.table(pheno,
            "GSE197674 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(beta,
            "GSE197674 beta after filtering.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

beta <- data.table::fread("GSE197674 beta after filtering.txt",
                          sep = "\t")
beta <- as.data.frame(beta)
rownames(beta) <- beta$V1
beta <- beta %>% select(-V1)
gc()
pheno <- read.delim("GSE197674 Phenotypes.txt")
pheno <- separate(pheno, col = supplementary_file, into = c("Suppinfo","Sentrix_ID", "Sentrix_Position","Colour"), sep = "\\_")
glimpse(pheno)

impute <- champ.impute(beta = as.matrix(beta),
                       pd = pheno)

#can't do champ.impute because do not have the processing power and R keeps crashing 

#row.has.na <- apply(beta, 1, function(x) any(is.na(x)))
#beta <- beta[!row.has.na,]

#this method removes too many rows 


#split the beta matrix into 2 and do imputations on wach subset 
beta1 <- beta[,c(1:1000)]
beta2 <- beta[,c(1001:2138)]
pheno1 <- pheno %>% filter(Sample_Name %in% colnames(beta1))
pheno2 <- pheno %>% filter(Sample_Name %in% colnames(beta2))
impute1 <- champ.impute(beta = as.matrix(beta1),
                        pd = pheno1)
impute2 <- champ.impute(beta = as.matrix(beta2),
                        pd = pheno2)
beta1 <- impute1$beta
beta1 <- as.data.frame(beta1)
beta2 <- impute2$beta
beta2 <- as.data.frame(beta2)

beta1$CpG <- rownames(beta1)
beta2$CpG <- rownames(beta2)
beta<- inner_join(beta1, beta2, by = "CpG")
rownames(beta) <- beta$CpG
CpG <- rownames(beta)
remove(beta1)
remove(beta2)
beta <- beta %>% select(-CpG)
rownames(beta)

write.table(pheno,
            "GSE197674 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(filter$beta,
            "GSE197674 beta after filtering and imputation.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

filter <- champ.filter(beta = as.matrix(beta),
                       pd = pheno, 
                       arraytype = "EPIC")

champ.QC(beta = filter$beta,
         pheno = filter$pd$sex,
         dendrogram = FALSE)

glimpse(pheno)
champ.SVD(beta=filter$beta,
          pd= pheno %>%select(Sentrix_ID,
                              Sentrix_Position,
                              age,
                              sex,
                              `abdominal_pelvic_rt.ch1`,
                              `alkylating.agent.classic.ch1`,
                              `anthracyclines.ch1`,
                              `brainrt.ch1`,
                              `chestrt.ch1`,
                              `corticosteroids.ch1`,
                              `epipodophyllotoxins.ch1`,
                              `platinum.ch1`,
                              `vincristine.ch1`),
          resultsDir="./CHAMP_SVDimages")

#try run combat to do a batch correction 

library(sva)
myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_ID,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

myCombat = logit2(myCombat)

myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=pheno$Sentrix_Position,
                mod=NULL)
gc()
myCombat=ilogit2(myCombat)

write.table(myCombat,
            "GSE197674 beta after filtering and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")



#### GSE197676 ####

path = "/pvol/Preprocessing/Blood/GEO/GSE197676"
setwd(path)

library(tidyverse)
library(ChAMP)
library(limma)
library(minfi)
library(data.table)
library(GEOquery)

#### Preprocessing 
getGEOSuppFiles("GSE197676",
                fetch_files = TRUE,
                filter_regex = "tar$")

pheno <- read.delim("GSE197676 Phenotypes.txt", sep = " ")
write.csv(pheno, "GSE197676 Phenotypes.csv")

untar("GSE197676/GSE197676_RAW.tar")

targets <- read.metharray.sheet(getwd())
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
detP <- detectionP(rgSet)
MSet <- preprocessRaw(rgSet) 
meth = getMeth(MSet)
unmeth = getUnmeth(MSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <- addSex(GRset, sex = predictedSex)

filter = champ.load(getwd(),
                    arraytype = "EPIC")

library(readxl)
sheets <- excel_sheets("/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = filter$beta
for (s in sheets)
{
  qcprobes <- read_excel('/pvol/Preprocessing/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}

pheno <- filter$pd
pheno <- as.data.frame(pheno)

champ.QC(beta = beta,
         pheno = pheno$sex,
         dendrogram = FALSE)

myNorm <- champ.norm(beta=beta,
                     arraytype="EPIC")

champ.QC(beta = myNorm,
         pheno = pheno$sex,
         dendrogram = FALSE)

glimpse(pheno)
champ.SVD(beta=myNorm,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex),
          resultsDir="./CHAMP_SVDimages/")

library(sva)
M <- logit2(myNorm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Slide),
                mod=NULL)

myCombat=ilogit2(myCombat)

glimpse(pheno)

champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex),
          resultsDir="./CHAMP_SVDimages/batchcorrected")

myCombat <- logit2(myCombat)
myCombat=ComBat(dat=as.matrix(myCombat), #it outputs an M-value matrix adjusted for batch
                batch=as.factor(pheno$Array),
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.QC(beta = myCombat,
         pheno = pheno$sex,
         dendrogram = FALSE)

champ.SVD(beta=myCombat,
          pd=pheno%>%select(Slide,
                            Array,
                            age,
                            sex),
          resultsDir="./CHAMP_SVDimages/batchandpositioncorrected")

write.table(pheno,
            "GSE197676 Phenotypes.txt",
            col.names = TRUE,
            sep = "\t")

write.table(myCombat,
            "GSE197676 beta after normalisation and batch correction.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

