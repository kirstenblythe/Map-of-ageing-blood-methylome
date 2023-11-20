#### Pathway enrichment ####

setwd("working directory")

library(tidyverse)
library(missMethyl)
library(tidyverse)

meta_total <- read.delim("Meta-analysis results DMPs blood.txt", sep = "")
colnames <- colnames(meta_total)
meta_total$Annotated_genes[is.na(meta_total$Annotated_genes)] <- "" #replace NA with nothing

#create our own annotation for each cpg 
RSanno <- DataFrame(chr = meta_total$Chromosome,
                    pos = meta_total$Position_hg38,
                    Name = meta_total$CpG,
                    UCSC_RefGene_Name = meta_total$Annotated_genes,
                    UCSC_RefGene_Group = meta_total$Annotated_genes,
                    row.names = meta_total$CpG)

#Run pathway analysis on MSigDB relevant gene sets
library(mitch)

#download genesets from https://www.gsea-msigdb.org/gsea/msigdb

genesets <- list.files(path = "/pvol/EWAS/Blood/MSigDB/msigdb_v2022.1.Hs_GMTs",
                       pattern = ".gmt",
                       recursive = TRUE)
cgp <- genesets[which(str_detect(genesets,
                                 fixed("cgp"))==TRUE)]
cp <- genesets[which(str_detect(genesets,
                                fixed("cp"))==TRUE)]
go <- genesets[which(str_detect(genesets,
                                fixed(".go."))==TRUE)]
hpo <- genesets[which(str_detect(genesets,
                                 fixed("hpo"))==TRUE)]
c7 <- genesets[which(str_detect(genesets,
                                fixed("c7"))==TRUE)]

setwd("/pvol/EWAS/Blood/MSigDB/msigdb_v2022.1.Hs_GMTs")
cgp <- gmt_import(cgp[which(str_detect(cgp,
                                       "entrez")==TRUE)])
cp <- gmt_import(cp[which(str_detect(cp,
                                     "biocarta|kegg|pid|reactome|wiki|symbols",
                                     negate = TRUE)==TRUE)])
go <- gmt_import(go[which(str_detect(go,
                                     ".bp|.cc|.mf|symbols",
                                     negate = TRUE)==TRUE)])
hpo <- gmt_import(hpo[which(str_detect(hpo,
                                       "entrez")==TRUE)])
c7 <- gmt_import(c7[which(str_detect(c7,
                                     "entrez")==TRUE&
                            str_detect(c7,
                                       "all")==TRUE)])

####Enrichment for each gene set

#First look at hyperDMPs vs hypoDMPs
DMPs <- read_csv("Blood_DMPs.csv") %>% dplyr::rename(CpG = MarkerName)
DMPs_hypo <- DMPs %>% filter(Effect < 0) 
DMPs_hyper <- DMPs %>% filter(Effect >= 0) 


#Now compare the three classes of age-related CpGs (homoDMPs, DMPs-VMPs, constantVMPs)
VMPs <- read_csv("/pvol/EWAS/Blood/Meta-analysis/VMPs/Blood_VMPs.csv") %>% dplyr::rename(CpG = MarkerName)

#Obtain the three classes of age-related CpGs
DMPs_only <- DMPs %>%
  filter(!CpG %in% VMPs$CpG)

DMPs_VMPs <- DMPs %>%
  filter(CpG %in% VMPs$CpG)

VMPs_only <- VMPs %>%
  filter(!CpG %in% DMPs$CpG)


#Remove gene sets < 10 genes
geneset <- cp
size <- sapply(geneset,length)
geneset <- geneset[-which(size<10|size>500)]


#######Is there a significant overlap of common pathways?

#set working directory where the msigDBs have been downloaded
setwd("/pvol/EWAS/Blood/MSigDB/msigdb_v2022.1.Hs_GMTs")
library(ActivePathways)
library(mitch)
library(tidyverse)
genesets <- list.files(path = "/pvol/EWAS/Blood/MSigDB/msigdb_v2022.1.Hs_GMTs",
                       pattern = ".gmt",
                       recursive = TRUE)
go <- genesets[which(str_detect(genesets,
                                fixed(".go."))==TRUE)]
cp <- genesets[which(str_detect(genesets,
                                fixed("cp"))==TRUE)]


go <- read.GMT(go[which(str_detect(go,
                                   ".bp|.cc|.mf|symbols",
                                   negate = TRUE)==TRUE)])
cp <- read.GMT(cp[which(str_detect(cp,
                                   "biocarta|kegg|pid|reactome|wiki|symbols",
                                   negate = TRUE)==TRUE)])

#enrichment of hyperDMPs vs hypoDMPs in each of the msig signatures beginning with CP enrichment
enrichment_hypo <- read_delim("CP enrichment hypoDMPs meta-analysis.txt") %>%
  select(ID,P.DE,FDR) %>%
  dplyr::rename(P.DE.hypo = P.DE,
                FDR.hypo = FDR)
enrichment_hyper <- read_delim("CP enrichment hyperDMPs meta-analysis.txt")%>%
  select(ID,P.DE,FDR) %>%
  dplyr::rename(P.DE.hyper = P.DE,
                FDR.hyper = FDR)
enrichment <- full_join(enrichment_hypo,
                        enrichment_hyper) %>%
  filter(FDR.hypo < 0.005 | FDR.hyper < 0.005) %>%
  filter(str_detect(ID,"REACT"))%>%
  dplyr::rename(term.id = ID) %>%
  mutate(adjusted.p.val = 0.001,
         term.name = term.id)
write.table(enrichment %>% select(term.id,term.name,adjusted.p.val),
            "enrichmentMap__pathways_REACT_hypohyper.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")
cp <- cp[enrichment %>% pull(term.id)]
write.GMT(cp,
          "enrichmentMap__pathways_REACT_hypohyper.gmt")

labels <- enrichment %>%
  mutate(hypo = ifelse(FDR.hypo<0.005,1,0),
         hyper = ifelse(FDR.hyper<0.005,1,0),
         instruct = 'piechart: attributelist="hypo,hyper" colorlist="#0000FF,#FF0000" showlabels=FALSE') %>%
  select(term.id,
         hypo,
         hyper,
         instruct)
write.table(labels,
            "enrichmentMap__subgroups_REAC_hypohyper.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Create Cytoscape files to create a network whose nodes are colored according to significance in hyper or hypo or both

#repeat for the homoDMPs, DMPs-VMPs and constant VMPs 
enrichment_homo <- read_delim("CP enrichment homoDMPs meta-analysis.txt") %>%
  select(ID,P.DE,FDR) %>%
  dplyr::rename(P.DE.homo = P.DE,
                FDR.homo = FDR)
enrichment_DMPVMP <- read_delim("CP enrichment DMPs-VMPs meta-analysis.txt") %>%
  select(ID,P.DE,FDR) %>%
  dplyr::rename(P.DE.DMPVMP = P.DE,
                FDR.DMPVMP = FDR)
enrichment_const <- read_delim("CP enrichment constant VMPs meta-analysis.txt")%>%
  select(ID,P.DE,FDR) %>%
  dplyr::rename(P.DE.const = P.DE,
                FDR.const = FDR)
enrichment <- Reduce(full_join,
                     list(enrichment_homo,
                          enrichment_DMPVMP,
                          enrichment_const)) %>%
  filter(FDR.homo < 0.005 | FDR.DMPVMP < 0.005 | FDR.const < 0.005) %>%
  filter(str_detect(ID,"REACT"))%>%
  dplyr::rename(term.id = ID) %>%
  mutate(adjusted.p.val = 0.001,
         term.name = term.id)
write.table(enrichment %>% select(term.id,term.name,adjusted.p.val),
            "enrichmentMap__pathways_REACT_DMPVMP.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#go <- go[enrichment %>% pull(term.id)]
cp <- cp[enrichment %>% pull(term.id)]
write.GMT(cp,
          "enrichmentMap__pathways_REACT_DMPVMP.gmt")

labels <- enrichment %>%
  mutate(homoDMP = ifelse(FDR.homo<0.005,1,0),
         DMPVMP = ifelse(FDR.DMPVMP<0.005,1,0),
         constVMP = ifelse(FDR.const<0.005,1,0),
         instruct = 'piechart: attributelist="homoDMP,DMPVMP,constVMP" colorlist="#f54242,#f542e3,#4287f5" showlabels=FALSE') %>%
  select(term.id,
         homoDMP,
         DMPVMP,
         constVMP,
         instruct)
write.table(labels,
            "enrichmentMap__subgroups_REAC_DMPVMP.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#####Chi-square test
matrix <- cbind(
  #Hypo
  c(nrow(enrichment %>% filter(FDR.hypo<0.005&FDR.hyper<0.005)),
    nrow(enrichment %>% filter(FDR.hypo<0.005&FDR.hyper>0.005))),
  #Hyper
  c(nrow(enrichment %>% filter(FDR.hypo>0.005&FDR.hyper<0.005)),
    nrow(enrichment %>% filter(FDR.hypo>0.005&FDR.hyper>0.005))))
colnames(matrix) <- c("Sig hypo",
                      "Non-sig hypo")
rownames(matrix) <- c("Sig hyper",
                      "Non-sig hyper")
chisq <- chisq.test(matrix)
chisq$stdres[1,]


enrichment_sig <- enrichment %>%
  filter(FDR.homo<0.005&FDR.DMPVMP<0.005&FDR.const<0.005)%>%
  filter(str_detect(ID,"GOBP"))
#filter(!str_detect(ID,"GSE")) 
#slice_head(n = 50)

#Compare lists of hypo vs hyper
sim <- term_similarity(gl = geneset[enrichment_sig$ID])
clust <- hclust(dist(sim), "average")
#clust$labels <- tolower(clust$labels)
clust$labels <- tolower(str_replace(clust$labels,
                                    "GOBP_",
                                    ""))

colors <- brewer.pal(n = 3, name = "Dark2")
clus4 = cutree(clust, 3)

tiff("GO gene set enrichment three classes of age-related CpGs.tiff",
     height = 5,
     width = 5,
     res = 300,
     unit = 'in')
plot(as.phylo(clust),
     type = "cladogram",
     tip.color = colors[clus4],
     cex = 0.6,
     no.margin = TRUE)
dev.off()

