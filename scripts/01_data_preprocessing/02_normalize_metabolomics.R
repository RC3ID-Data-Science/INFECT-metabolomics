## Install Metaboanalyst R#################

metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

BiocManager::install("ctc")
BiocManager::install("glasso")
BiocManager::install("huge")
BiocManager::install("ppcor")

install.packages(pacman)
library(pacman)
pacman::p_load(tidyverse,ctc,glasso, huge, ppcor,plotly,impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, SSPA, sva, limma, KEGGgraph, siggenes,BiocParallel, MSnbase, multtest,RBGL,edgeR,fgsea,httr,qs)

install.packages("usethis")
library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
#if still cannot be installed, download the file here https://drive.google.com/file/d/1OHRUzXFDukWXEKxTLMF9X1fHRN06uFFu/view
# then, go to terminal (beside console) and write...
#cd ~/Downloads
#R CMD INSTALL MetaboAnalystR_3.0.3.tar.gz


## load Metaboanalyst R #############
pacman::p_load(tidyverse,scales,dplyr,ctc,glasso, huge, ppcor,plotly,impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, SSPA, sva, limma, KEGGgraph, siggenes,BiocParallel, MSnbase, multtest,RBGL,edgeR,fgsea,httr,qs)
library(MetaboAnalystR)

sim_HM_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/2
    x
  })
  return(result)
}

## file name and location #########
dir     <- "D:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\"

HILIC_pos   <- paste0(dir,"1. HILIC-pos_annotated.HMDB.csv")
HILIC_neg   <- paste0(dir,"1. HILIC-neg_annotated.HMDB.csv")
C8      <- paste0(dir,"1. C8-pos_annotated.HMDB.csv")
C18     <- paste0(dir,"1. C18-neg_annotated.HMDB.csv")

dataset<- c(HILIC_pos=HILIC_pos, 
            HILIC_neg = HILIC_neg, 
            C8=C8, 
            C18=C18)

# 1. Check distribution of missing value from all dataset ---------------------------------

list2env(sapply(dataset,read.csv),envir=.GlobalEnv)

missing <- NULL

for (y in dataset) {
  x <- read.csv(y)
  status <- x[1:2]
  x <- x %>% 
   group_by(status) %>%
    select(everything()) %>%
    summarise_all(funs(sum(is.na(.)))) %>% as_tibble()
  x <- as.data.frame(t(x))
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1:2),]
  x$Converter <- as.numeric(x$Converter)
  x$Persistently_uninfected <- as.numeric(x$Persistently_uninfected)
  x$total <- x$Converter+x$Persistently_uninfected
  x<-filter(x,total!=0)
  x$per.converter =percent(x$Converter/71)
  x$per.persistent=percent(x$Persistently_uninfected/128)
  x$per.total=percent((x$total)/199)
  
  missing <- as.data.frame(rbind(missing,x))
}

status$samples <- as.character(status$samples)
status <- status %>% dplyr::rename(IdInfect = samples)

#view 'missing' to look at the distribution of missing value. think the cut-off point by looking at the total missing percentage for omitting variables in the next step#

# 2. hydrophilic interaction liquid chromatography (HILIC)-positive (AAs, AA metabolites, acylcarnitines, dipeptides, and other cationic polar metabolites) data processing ------------------------------------------------
# * 2.1. Delete variables and replace missing value ------------------------
HILIC_pos   <- paste0(dir,"1. HILIC-pos_annotated.HMDB.csv")



mSet_HILIC_pos <- InitDataObjects("list", "stat", FALSE)
mSet_HILIC_pos <- Read.TextData(mSet_HILIC_pos,HILIC_pos, "row", "disc")
mSet_HILIC_pos <- SanityCheckData(mSet_HILIC_pos)
mSet_HILIC_pos <- RemoveMissingPercent(mSetObj = mSet_HILIC_pos,percent = 0.15)
# 15% is choosen after looking at the distribution of missing value. metabolites that have more than 15% missing value will be omitted#
mSet_HILIC_pos$msgSet$replace.msg
mSet_HILIC_pos$dataSet$proc <- sim_HM_wrapper(mSet_HILIC_pos$dataSet$preproc)
#mSet_HILIC_pos <- ReplaceMin(mSet_HILIC_pos)
#mSet_HILIC_pos <- ImputeVar(mSet_HILIC_pos,method = "min")
mSet_HILIC_pos$msgSet$replace.msg

# * 2.2. normalization and transformation ----------------------------------

#mSet_HILIC <- FilterVariable(mSetObj = mSet_HILIC_pos,filter = "iqr",qcFilter = F)
#?FilterVariable
mSet_HILIC_pos<-PreparePrenormData(mSet_HILIC_pos)
colnames(mSet_HILIC_pos$dataSet$prenorm)
#feature.nm.vec <- c("")
#smpl.nm.vec <- c("")
#grp.nm.vec <- c("QC")
#mSet <- UpdateData(mSet)
#mSet_HILIC<-Normalization(mSet_HILIC, rowNorm = "None", transNorm = "None", scaleNorm = "None")
mSet_HILIC_pos <-Normalization(mSet_HILIC_pos, rowNorm = "None", transNorm = "LogNorm", scaleNorm = "None")
mSet_HILIC_pos <- PlotNormSummary(mSet_HILIC_pos, "norm_HILIC_pos", "png", 72, width=NA)
mSet_HILIC_pos <- PlotSampleNormSummary(mSet_HILIC_pos, "norm_HILIC_pos.Sample_", "png", 72, width=NA)
# * 2.3. HILIC processed datasets -----------------------------------------------------

HILIC_pos.data <- list(normalized=mSet_HILIC_pos$dataSet$norm,processed=mSet_HILIC_pos$dataSet$proc,original=mSet_HILIC_pos$dataSet$orig)

#finding the omitted variables
ori_meta <- c(colnames(HILIC_pos.data$original))
proc_meta <- c(colnames(HILIC_pos.data$processed))
norm_meta <- c(colnames(HILIC_pos.data$normalized))
g <- as.data.frame(ori_meta)
omitted.proc.meta.HILIC <- g %>% filter(!g$ori_meta %in% proc_meta)
omitted.norm.meta.HILIC <- g %>% filter(!g$ori_meta %in% norm_meta)


write.csv(HILIC_pos.data$normalized,paste0(dir,"1. HILIC-pos_normalized.HMDB.csv"),row.names = FALSE)
write.csv(HILIC_pos.data$processed,paste0(dir,"1. HILIC-pos_processed.HMDB.csv"),row.names = FALSE) 
write.csv(HILIC_pos.data$original,paste0(dir,"1. HILIC-pos_raw.HMDB.csv"),row.names = FALSE)

# 2.b HILIC neg -----------------------------------------------------------

# 2. hydrophilic interaction liquid chromatography (HILIC)-positive (AAs, AA metabolites, acylcarnitines, dipeptides, and other cationic polar metabolites) data processing ------------------------------------------------
# * 2.1. Delete variables and replace missing value ------------------------


HILIC_neg   <- paste0(dir,"1. HILIC-neg_annotated.HMDB.csv")



mSet_HILIC_neg <- InitDataObjects("list", "stat", FALSE)
mSet_HILIC_neg <- Read.TextData(mSet_HILIC_neg,HILIC_neg, "row", "disc")
mSet_HILIC_neg <- SanityCheckData(mSet_HILIC_neg)
mSet_HILIC_neg <- RemoveMissingPercent(mSetObj = mSet_HILIC_neg,percent = 0.15)
# 15% is choosen after looking at the distribution of missing value. metabolites that have more than 15% missing value will be omitted#
mSet_HILIC_neg$msgSet$replace.msg
mSet_HILIC_neg$dataSet$proc <- sim_HM_wrapper(mSet_HILIC_neg$dataSet$preproc)
#mSet_HILIC_neg <- ReplaceMin(mSet_HILIC_neg)
#mSet_HILIC_neg <- ImputeVar(mSet_HILIC_neg,method = "min")
mSet_HILIC_neg$msgSet$replace.msg

# * 2.2. normalization and transformation ----------------------------------

#mSet_HILIC <- FilterVariable(mSetObj = mSet_HILIC_neg,filter = "iqr",qcFilter = F)
#?FilterVariable
mSet_HILIC_neg<-PreparePrenormData(mSet_HILIC_neg)
colnames(mSet_HILIC_neg$dataSet$prenorm)
#feature.nm.vec <- c("")
#smpl.nm.vec <- c("")
#grp.nm.vec <- c("QC")
#mSet <- UpdateData(mSet)
#mSet_HILIC<-Normalization(mSet_HILIC, rowNorm = "None", transNorm = "None", scaleNorm = "None")
mSet_HILIC_neg <-Normalization(mSet_HILIC_neg, rowNorm = "None", transNorm = "LogNorm", scaleNorm = "None")
mSet_HILIC_neg <- PlotNormSummary(mSet_HILIC_neg, "norm_HILIC_neg", "png", 72, width=NA)
mSet_HILIC_neg <- PlotSampleNormSummary(mSet_HILIC_neg, "norm_HILIC_neg.Sample_", "png", 72, width=NA)
# * 2.3. HILIC processed datasets -----------------------------------------------------

HILIC_neg.data <- list(normalized=mSet_HILIC_neg$dataSet$norm,processed=mSet_HILIC_neg$dataSet$proc,original=mSet_HILIC_neg$dataSet$orig)

#finding the omitted variables
ori_meta <- c(colnames(HILIC_neg.data$original))
proc_meta <- c(colnames(HILIC_neg.data$processed))
norm_meta <- c(colnames(HILIC_neg.data$normalized))
g <- as.data.frame(ori_meta)
omitted.proc.meta.HILIC <- g %>% filter(!g$ori_meta %in% proc_meta)
omitted.norm.meta.HILIC <- g %>% filter(!g$ori_meta %in% norm_meta)


write.csv(HILIC_neg.data$normalized,paste0(dir,"1. HILIC-neg_normalized.HMDB.csv"),row.names = FALSE)
write.csv(HILIC_neg.data$processed,paste0(dir,"1. HILIC-neg_processed.HMDB.csv"),row.names = FALSE) 
write.csv(HILIC_neg.data$original,paste0(dir,"1. HILIC-neg_raw.HMDB.csv"),row.names = FALSE)
# 3. C18-negative (free fatty acids, eicosanoids, bile acids, and metabolites of intermediate polarity) data processing ------------------------------------------------
# * 3.1. Delete variables and replace missing value ------------------------


C18     <- paste0(dir,"1. C18-neg_annotated.HMDB.csv")

mSet_C18 <- InitDataObjects("list", "stat", FALSE)
mSet_C18 <- Read.TextData(mSet_C18,C18, "row", "disc")
mSet_C18 <- SanityCheckData(mSet_C18)
mSet_C18 <- RemoveMissingPercent(mSetObj = mSet_C18,percent = 0.15)
mSet_C18$msgSet$replace.msg
mSet_C18$dataSet$proc <- sim_HM_wrapper(mSet_C18$dataSet$preproc)
#mSet_C18 <- ReplaceMin(mSet_C18)
#mSet_C18 <- ImputeVar(mSet_C18,method = "min")
mSet_C18$msgSet$replace.msg

# * 3.2. normalization and transformation ----------------------------------

#mSet_C18 <- FilterVariable(mSetObj = mSet_C18,filter = "iqr",qcFilter = F)
#?FilterVariable
mSet_C18<-PreparePrenormData(mSet_C18)
colnames(mSet_C18$dataSet$prenorm)
#feature.nm.vec <- c("")
#smpl.nm.vec <- c("")
#grp.nm.vec <- c("QC")
#mSet <- UpdateData(mSet)
mSet_C18<-Normalization(mSet_C18, rowNorm = "None", transNorm = "LogNorm", scaleNorm = "None")
#mSet_C18<-Normalization(mSet_C18, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "15R-15-methyl PGA2_IS")
#mSet_C18<-Normalization(mSet_C18, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "15S-15-methyl PGD2 iSTD")
#mSet_C18<-Normalization(mSet_C18, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "15S-15-methyl PGE1_IS")
#mSet_C18<-Normalization(mSet_C18, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "15S-15-methyl PGE2_IS")
#mSet_C18<-Normalization(mSet_C18, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "15R-15-methyl PGF2alpha_IS")
mSet_C18 <- PlotNormSummary(mSet_C18, "norm_C18_", "png", 72, width=NA)
mSet_C18 <- PlotSampleNormSummary(mSet_C18, "mSet_C18.Sample_", "png", 72, width=NA)

# * 3.3. C18 processed datasets -----------------------------------------------------

C18.data <- list(normalized=mSet_C18$dataSet$norm,processed=mSet_C18$dataSet$proc,original=mSet_C18$dataSet$orig)

#finding the omitted variables
ori_meta <- c(colnames(C18.data$original))
proc_meta <- c(colnames(C18.data$processed))
norm_meta <- c(colnames(C18.data$normalized))
g <- as.data.frame(ori_meta)
omitted.proc.meta.C18 <- g %>% filter(!g$ori_meta %in% proc_meta)
omitted.norm.meta.C18 <- g %>% filter(!g$ori_meta %in% norm_meta)

write.csv(C18.data$normalized,paste0(dir,"1. C18-neg_normalized.HMDB.csv"), row.names = FALSE)
write.csv(C18.data$processed,paste0(dir,"1. C18-neg_processed.HMDB.csv"), row.names = FALSE)
write.csv(C18.data$original,paste0(dir,"1. C18-neg_raw.HMDB.csv"), row.names = FALSE)


# 4. C8-positive (lipids and nonpolar metabolites lipids) data processing ------------------------------------------------
# * 4.1. Delete variables and replace missing value ------------------------


C8      <- paste0(dir,"1. C8-pos_annotated.HMDB.csv")


mSet_C8 <- InitDataObjects("list", "stat", FALSE)
mSet_C8 <- Read.TextData(mSet_C8,C8, "row", "disc")
mSet_C8 <- SanityCheckData(mSet_C8)
mSet_C8 <- RemoveMissingPercent(mSetObj = mSet_C8,percent = 0.15)
mSet_C8$msgSet$replace.msg
mSet_C8$dataSet$proc <- sim_HM_wrapper(mSet_C8$dataSet$preproc)
#mSet_C8 <- ReplaceMin(mSet_C8)
#mSet_C8 <- ImputeVar(mSet_C8,method = "min")
mSet_C8$msgSet$replace.msg

# * 4.2. normalization and transformation ----------------------------------

#mSet_C8 <- FilterVariable(mSetObj = mSet_C8,filter = "iqr",qcFilter = F)
#?FilterVariable
mSet_C8<-PreparePrenormData(mSet_C8)
colnames(mSet_C8$dataSet$prenorm)
#feature.nm.vec <- c("")
#smpl.nm.vec <- c("")
#grp.nm.vec <- c("QC")
#mSet <- UpdateData(mSet)
#mSet_C8<-Normalization(mSet_C8, rowNorm = "CompNorm", transNorm = "LogNorm", scaleNorm = "None", ref = "PC120/120")
mSet_C8 <-Normalization(mSet_C8, rowNorm = "None", transNorm = "LogNorm", scaleNorm = "None")
mSet_C8 <- PlotNormSummary(mSet_C8, "norm_C8_", "png", 72, width=NA)
mSet_C8 <- PlotSampleNormSummary(mSet_C8, "mSet_C8.Sample_", "png", 72, width=NA)

# * 4.3. C8 processed datasets -----------------------------------------------------

C8.data <- list(normalized=mSet_C8$dataSet$norm,processed=mSet_C8$dataSet$proc,original=mSet_C8$dataSet$orig)

#finding the omitted variables
ori_meta <- c(colnames(C8.data$original))
proc_meta <- c(colnames(C8.data$processed))
norm_meta <- c(colnames(C8.data$normalized))
g <- as.data.frame(ori_meta)
omitted.proc.meta.C8 <- g %>% filter(!g$ori_meta %in% proc_meta)
omitted.norm.meta.C8 <- g %>% filter(!g$ori_meta %in% norm_meta)


write.csv(C8.data$normalized,paste0(dir,"1. C8-pos_normalized.HMDB.csv"), row.names = FALSE) 
write.csv(C8.data$processed,paste0(dir,"1. C8-pos_processed.HMDB.csv"), row.names = FALSE) 
write.csv(C8.data$original,paste0(dir,"1. C8-pos_raw.HMDB.csv"), row.names = FALSE)


# 5. Merge all normalized datasets ----------------------------------------

HILIC_pos_norm <- as.data.frame(HILIC_pos.data$normalized)
HILIC_pos_norm <- tibble::rownames_to_column(HILIC_pos_norm, "IdInfect")
HILIC_neg_norm <- as.data.frame(HILIC_neg.data$normalized)
HILIC_neg_norm <- tibble::rownames_to_column(HILIC_neg_norm, "IdInfect")
C18_norm <- as.data.frame(C18.data$normalized)
C18_norm <- tibble::rownames_to_column(C18_norm, "IdInfect")
C8_norm <- as.data.frame(C8.data$normalized)
C8_norm <- tibble::rownames_to_column(C8_norm, "IdInfect")


all_norm_meta <- full_join(HILIC_pos_norm,HILIC_neg_norm, by="IdInfect", suffix = c(".HIL_pos",".HIL_neg"))
all_norm_meta <- full_join(all_norm_meta,C8_norm, by="IdInfect", suffix = c("",".C8"))
all_norm_meta <- full_join(all_norm_meta,C18_norm, by="IdInfect", suffix = c("",".C18"))
all_norm_meta <- full_join(all_norm_meta,status,by = "IdInfect")

#all_norm_meta <- merge(HILIC_pos,HILIC_neg, by = "IdInfect", all.x = TRUE)
#all_norm_meta <- merge(all_norm_meta,C8_norm, by = "IdInfect", all.x = TRUE)  
#all_norm_meta <- merge(all_norm_meta,status,by = "IdInfect", all.x=TRUE)
all_norm_meta <- all_norm_meta %>% relocate(status, .after = IdInfect)
all_norm_meta <- all_norm_meta %>% dplyr::rename( samples = IdInfect)

all_norm_meta <- all_norm_meta %>% select(!starts_with(c("internal","NA","redundant")))

getCV <- function(df, minimum_cv=0, ...){
  
  
  cv_qc <- apply(df, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100})
  cv_qc <- cv_qc[cv_qc >= minimum_cv]
  ord <- order(cv_qc)
  
  barplot(cv_qc[ord], 
          horiz = T, las=1, xlim = c(0,100), ...)
  abline(v=minimum_cv, lty=2, col='red')
  
  return(names(cv_qc))
}

high_cv_feature <- getCV(all_norm_meta[,-c(1:2)], minimum_cv = 30, cex.names = 0.6)
high_cv_feature

cv_qc <- apply(all_norm_meta[,-c(1:2)], 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100})
cv_qc

hist(cv_qc)


all_norm_meta <- all_norm_meta[,!(colnames(all_norm_meta) %in% high_cv_feature)]



write.csv(all_norm_meta,paste0(dir,"1.all_norm_metabolomics.HMDB.csv"), row.names = FALSE)

###################################------------------------

cor(all_norm_meta$HMDB0012103,all_norm_meta$HMDB0012103.C8)


all_norm_meta_analysis <- all_norm_meta[,-2]
write.csv(all_norm_meta_analysis,paste0(dir,"1.all_norm_metabolomics_only.HMDB.csv"),row.names = FALSE)

dir     <- "C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\"
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, paste0(dir,"1.all_norm_metabolomics.HMDB.csv"), "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, pcx = 1, pcy = 2, reg = 0.95, show=0 ,grey.scale = 0)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_sample_", "png", 72, width=NA, pcx = 1, pcy = 2, reg = 0.95, show=1 ,grey.scale = 0)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)


dir     <- "C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\"
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet,paste0(dir,"1. all_norm_metabolomics.csv"), "rowu", "disc")
mSet <- SanityCheckData(mSet)
#mSet <- PreparePrenormData(mSet)
#feature.nm.vec <- c("")
#smpl.nm.vec <- c("")
#grp.nm.vec <- c("QC")
#mSet <- UpdateData(mSet)
#mSet <- Normalization(mSet, rowNorm = "None", transNorm = "LogNorm", scaleNorm = "None")
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_HILIC_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_HILIC_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_HILIC_", "png", 72, width=NA, pcx = 1, pcy = 2, reg = 0.95, 0,0)

mSet<-FC.Anal.unpaired(mSetObj = mSet,fc.thresh = 1.2, cmp.type = 1)
mSet <- PlotFC(mSetObj = mSet,"FC_HILIC_","png",72,width = NA)
mSet<-Ttests.Anal(mSetObj = mSet,nonpar = F, threshp = 0.05, paired = FALSE, equal.var = TRUE, all_results = TRUE)

?Ttests.Anal

mSet<-PlotTT(mSet, "tt1_HILIC_", "png", 72, width=NA)

z <- data.frame(mSet$analSet$tt$p.value)

PlotCmpdView(mSet, "Glutamic acid", "png", 72, width=NA)
PlotCmpdView(mSet, "Serotonin", "png", 72, width=NA)
PlotCmpdView(mSet, "Proline", "png", 72, width=NA)
PlotCmpdView(mSet, "N-Acetylaspartic acid", "png", 72, width=NA)
PlotCmpdView(mSet, "Trimethylamine N-oxide", "png", 72, width=NA)
PlotCmpdView(mSet, "Tryptophan", "png", 72, width=NA)
PlotCmpdView(mSet, "PSP-361/PSO-362", "png", 72, width=NA)


mSet<-Volcano.Anal(mSetObj = mSet, paired = FALSE, fcthresh = 1.2, cmpType = 1, F, threshp =  0.05, equal.var= TRUE, pval.type = "raw")
mSet<-PlotVolcano(mSet, "volcano_HILIC_", 1, "png", 300, width=NA)
mSet<-Volcano.Anal(mSetObj = mSet, paired = FALSE, fcthresh = 1.2, cmpType = 1, F, threshp =  0.5, TRUE, pval.type = "fdr")
mSet<-PlotVolcano(mSet, "volcano_HILIC_fdr",1, "png", 300, width=7.2)

?PlotVolcano
?Volcano.Anal
?PlotCmpdSummary

mSet<-PlotCorrHeatMap(mSetObj=mSet, "corr_overview_", "png", 300, width=NA, target="col", cor.method =  "pearson", colors =  "bwm", viewOpt = "overview", fix.col =  FALSE, no.clst =  FALSE, top = FALSE, topNum = "0")
mSet<-PlotCorrHeatMap(mSetObj=mSet, "corr_detail_", "png", 300, width=NA, target="col", cor.method =  "pearson", colors =  "bwm", viewOpt = "detail", fix.col =  FALSE, no.clst =  FALSE, top = FALSE, topNum = "0")
x<-read.csv("C:\\Users\\todia\\Google Drive\\1. PhD\\INFECT\\Metabolomics\\working\\HILIC_3_BC.csv")

str(x)

x$samples<-sapply(x$samples,as.character)

z<-raw %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.)))) %>% as.tibble()

raw

write.csv(z,"C:\\Users\\todia\\Google Drive\\1. PhD\\INFECT\\Metabolomics\\working\\HILIC_missing.csv")              

a <- as.data.frame(t(z))
za<-data.frame(colSums(a[,3:244]))
colsums
a<-filter(a,V1!=0)

write.csv(zb,"C:\\Users\\todia\\Google Drive\\1. PhD\\INFECT\\Metabolomics\\working\\HILIC_missing2.csv") 

a<-read.csv("C:\\Users\\todia\\Google Drive\\1. PhD\\INFECT\\Metabolomics\\working\\HILIC_3_BC.csv")


nama<-rownames(orig)
zc <- orig %>% 
  select(c(nama))
zc <- data.frame(t(zc))
zc <-`colnames<-`(zc,zc[1,])
zc <- zc[-c(1),]
str(zc)
zc<-filter(zc,)
zc[1:2]<-sapply(zc[1:2],as.numeric)
zc$per.converter =percent(zc$Converter/71)
zc$per.persistent=percent(zc$`Persistently uninfected`/128)

a$per.total=percent((a$V1)/199)
zc$omit = zc$per.total>percent(40)

zc <- filter(zc,per.total>50)
write.csv(a,"C:\\Users\\todia\\Google Drive\\1. PhD\\INFECT\\Metabolomics\\working\\HILIC_missing3.csv") 

colnames(mSet$dataSet$preproc)

colnames(mSet$dataSet$orig)
