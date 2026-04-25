# # load packages ---------------------------------------------------------
library(parallel)
library(pacman)
pacman::p_load(tidyverse, readxl, psych, devtools, openxlsx, limma, ggfortify, janitor)
pacman::p_load(dplyr, devtools,openxlsx,data.table,tableone,tidyverse)
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(pacman)
  library(dplyr)
  library(tidyverse)
  library(readxl)
})

# # File name and locations -----------------------------------------------
dir <- "D:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\"

# (1) metabolomics dataset
HIL_pos     <- "1. 20_0401_INFECT_HILIC-pos_metabolomics_results.xlsx"
HIL_neg     <- "1. 21_0305_INFECT_HILIC-neg_metabolomics_results.xlsx"
C8_pos      <- "1. 20_0401_INFECT_C8-pos_metabolomics_results.xlsx"
C18_neg     <- "1. 20_0401_INFECT_C18-neg_metabolomics_results.xlsx"
outcome     <- "2. assigned outcome.xlsx"

# # Open dataset ------------------------------------------------------------
list2env(parSapply(cl,list(HIL_pos = paste0(dir,HIL_pos), 
                           HIL_neg = paste0(dir,HIL_neg),
                           C8_pos  = paste0(dir,C8_pos), 
                           C18_neg = paste0(dir,C18_neg), 
                           outcome = paste0(dir,outcome)),
                read_excel),envir=.GlobalEnv)

{outcome <- outcome[,-c(2,4)]
outcome <- outcome %>% dplyr::rename(samples = IdInfect, status = outcome2)
outcome$samples <- as.character(outcome$samples)
clusterExport(cl, "outcome")
}

# 1. Making compound ID as variables. useful for analysis of unannotated metabolites --------

CompoundID <- list(HIL_pos = HIL_pos, HIL_neg = HIL_neg, C8_pos = C8_pos,C18_neg = C18_neg)
com <- function(x) {
  x <- x[-c(1:4),-c(1,3:7)]
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- as.data.frame(t(x))
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- tibble::rownames_to_column(x, "samples")
  x <- merge(x, outcome, by = "samples", all.x = TRUE)
  x <- x %>% relocate(status, .after = samples)
  x <- x[complete.cases(x[,2]), ]
  x <- x %>% 
  mutate(status = replace(status, status == 1, "Converter"))
  x <- x %>% 
  mutate(status = replace(status, status == 2, "Persistently_uninfected"))
}
CompoundID <- parSapply(cl,CompoundID, com)  

HIL_pos_MZ.RT <- HIL_pos[-c(1:4),c(1:4,7)]
HIL_pos_MZ.RT <- HIL_pos_MZ.RT %>% row_to_names(row_number = 1)
HIL_pos_MZ.RT$Compound_ID <- HIL_pos_MZ.RT$Compound_ID %>% paste0("_HILIC_pos")
HIL_pos_MZ.RT$Method <- "positive"
HIL_pos_MZ.RT <- HIL_pos_MZ.RT %>%
  mutate(Compound  = ifelse(is.na(HIL_pos_MZ.RT$Metabolite), HIL_pos_MZ.RT$Compound_ID, HIL_pos_MZ.RT$Metabolite))
HIL_pos_MZ.RT$Metabolite <- NULL


HIL_neg_MZ.RT <- HIL_neg[-c(1:4),c(1:4,7)]
HIL_neg_MZ.RT <- HIL_neg_MZ.RT %>% row_to_names(row_number = 1)
HIL_neg_MZ.RT$Compound_ID <- HIL_neg_MZ.RT$Compound_ID %>% paste0("_HILIC_neg")
HIL_neg_MZ.RT$Method <- "negative"
HIL_neg_MZ.RT <- HIL_neg_MZ.RT %>%
  mutate(Compound  = ifelse(is.na(HIL_neg_MZ.RT$Metabolite), HIL_neg_MZ.RT$Compound_ID, HIL_neg_MZ.RT$Metabolite))
HIL_neg_MZ.RT$Metabolite <- NULL


C8_pos_MZ.RT <- C8_pos[-c(1:4),c(1:4,7)]
C8_pos_MZ.RT <- C8_pos_MZ.RT %>% row_to_names(row_number = 1)
C8_pos_MZ.RT$Compound_ID <- C8_pos_MZ.RT$Compound_ID %>% paste0("_C8")
C8_pos_MZ.RT$Method <- "positive"
C8_pos_MZ.RT <- C8_pos_MZ.RT %>%
  mutate(Compound  = ifelse(is.na(C8_pos_MZ.RT$Metabolite), C8_pos_MZ.RT$Compound_ID, C8_pos_MZ.RT$Metabolite))
C8_pos_MZ.RT$Metabolite <- NULL

C18_neg_MZ.RT <- C18_neg[-c(1:4),c(1:4,7)]
C18_neg_MZ.RT <- C18_neg_MZ.RT %>% row_to_names(row_number = 1)
C18_neg_MZ.RT$Compound_ID <- C18_neg_MZ.RT$Compound_ID %>% paste0("_C18")
C18_neg_MZ.RT$Method <- "negative"
C18_neg_MZ.RT <- C18_neg_MZ.RT %>%
  mutate(Compound  = ifelse(is.na(C18_neg_MZ.RT$Metabolite), C18_neg_MZ.RT$Compound_ID, C18_neg_MZ.RT$Metabolite))
C18_neg_MZ.RT$Metabolite <- NULL


write.csv(HIL_pos_MZ.RT,paste0(dir,"1. HIL_pos_MZ.RT.csv"))
write.csv(HIL_neg_MZ.RT,paste0(dir,"1. HIL_neg_MZ.RT.csv"))
write.csv(C8_pos_MZ.RT,paste0(dir,"1. C8_pos_MZ.RT.csv"))
write.csv(C18_neg_MZ.RT,paste0(dir,"1. C18_neg_MZ.RT.csv"))

# change to factor and numeric (can be optimized)
  CompoundID$HIL_pos$status  <- as.factor(CompoundID$HIL_pos$status)
  CompoundID$HIL_neg$status  <- as.factor(CompoundID$HIL_neg$status)
  CompoundID$C8_pos$status   <- as.factor(CompoundID$C8_pos$status)
  CompoundID$C18_neg$status  <- as.factor(CompoundID$C18_neg$status)
  CompoundID$HIL_pos[,-c(1:2)] <- parSapply(cl,CompoundID$HIL_pos[,-c(1:2)],as.numeric)
  CompoundID$HIL_neg[,-c(1:2)] <- parSapply(cl,CompoundID$HIL_neg[,-c(1:2)],as.numeric)
  CompoundID$C8_pos[,-c(1:2)]  <- parSapply(cl,CompoundID$C8_pos[,-c(1:2)],as.numeric)
  CompoundID$C18_neg[,-c(1:2)] <- parSapply(cl,CompoundID$C18_neg[,-c(1:2)],as.numeric)

  
# # 1.1. save compound ID for unannotated -------------------------------------------------

write.csv(CompoundID$HIL_pos,paste0(dir,"1. HILIC-pos_X.unannotated.csv"), row.names = FALSE)
write.csv(CompoundID$HIL_neg,paste0(dir,"1. HILIC-neg_X.unannotated.csv"), row.names = FALSE)
write.csv(CompoundID$C8_pos,paste0(dir,"1. C8-pos_X.unannotated.csv"), row.names = FALSE)
write.csv(CompoundID$C18_neg,paste0(dir,"1. C18-neg_X.unannotated.csv"), row.names = FALSE)

# 2. Making metabolites names as the variables. useful for analysis --------
#reopening datasets

Metabolites <- list(HIL_pos = HIL_pos, HIL_neg = HIL_neg, C8_pos = C8_pos,C18_neg = C18_neg)


met <- function (x) {
x <- x[-c(1:4),-c(1:6)]
x <-`colnames<-`(x,x[1,])
x <- x[-c(1),]
x <- x[complete.cases(x[,1]), ]
x[,-1] <- sapply(x[,-1],as.numeric)
x <- as.data.frame(t(x))
x <-`colnames<-`(x,x[1,])
x <- x[-c(1),]
x <- tibble::rownames_to_column(x, "samples")
x <- merge(x, outcome, by = "samples", all.x = TRUE)
x <- x %>% relocate(status, .after = samples)
x <- x[complete.cases(x[,2]), ]
x <- x %>% 
  mutate(status = replace(status, status == 1, "Converter"))
x <- x %>% 
  mutate(status = replace(status, status == 2, "Persistently_uninfected"))
}
  
Metabolites <- parSapply(cl,Metabolites,met)



Metabolites$HIL_pos$status  <- as.factor(Metabolites$HIL_pos$status)
Metabolites$HIL_neg$status  <- as.factor(Metabolites$HIL_neg$status)
Metabolites$C8_pos$status   <- as.factor(Metabolites$C8_pos$status)
Metabolites$C18_neg$status  <- as.factor(Metabolites$C18_neg$status)
Metabolites$HIL_pos[,-c(1:2)] <- parSapply(cl,Metabolites$HIL_pos[,-c(1:2)],as.numeric)
Metabolites$HIL_neg[,-c(1:2)] <- parSapply(cl,Metabolites$HIL_neg[,-c(1:2)],as.numeric)
Metabolites$C8_pos[,-c(1:2)]  <- parSapply(cl,Metabolites$C8_pos[,-c(1:2)],as.numeric)
Metabolites$C18_neg[,-c(1:2)] <- parSapply(cl,Metabolites$C18_neg[,-c(1:2)],as.numeric)


# # 2.1. save metabolites names for annotated -------------------------------------------------

write.csv(Metabolites$HIL_pos,paste0(dir,"1. HILIC-pos_annotated.csv"), row.names = FALSE)
write.csv(Metabolites$HIL_neg,paste0(dir,"1. HILIC-neg_annotated.csv"), row.names = FALSE)
write.csv(Metabolites$C8_pos,paste0(dir,"1. C8-pos_annotated.csv"), row.names = FALSE)
write.csv(Metabolites$C18_neg,paste0(dir,"1. C18-neg_annotated.csv"), row.names = FALSE)


#  3. Making HMDB as the variables---------------------

HMDB_ID <- list(HIL_pos = HIL_pos, HIL_neg = HIL_neg, C8_pos = C8_pos,C18_neg = C18_neg)

hmdb_compound <- function(x) {
  x <- x[-c(1:4),-c(1:4)]
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- x[complete.cases(x[,1]), ]
  x <- x %>% select(c("HMDB_ID","Metabolite"))
  x <- x %>% filter(!str_detect(HMDB_ID,"NA|redundant|internal|REMOVE|Internal"))
}



HIL_pos <- hmdb_compound(HIL_pos)
HIL_pos$HMDB_ID <- gsub("\\*","",HIL_pos$HMDB_ID)
HIL_neg <- hmdb_compound(HIL_neg)
HIL_neg$HMDB_ID <- gsub("\\*","",HIL_neg$HMDB_ID)
C18_neg <- hmdb_compound(C18_neg)
C18_neg$HMDB_ID <- gsub("\\*","",C18_neg$HMDB_ID)
C8_pos  <- hmdb_compound(C8_pos)
C8_pos$HMDB_ID <- gsub("\\*","",C8_pos$HMDB_ID)

xx.1 <- list(HIL_pos = HIL_pos$HMDB_ID, HIL_neg = HIL_neg$HMDB_ID, C18_neg = C18_neg$HMDB_ID, C8_pos = C8_pos$HMDB_ID)

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

combs <- 
  unlist(lapply(1:length(xx.1), 
                function(j) combn(names(xx.1), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)
combs

elements <- 
  lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

Setdiff(xx.1[c("C18_neg", "C8_pos")], xx.1[c("HIL_neg", "HIL_pos")])
#[1] "O" "P" "K" "H"
#("HIL_pos", "C8_pos")], xx.1[c("HIL_neg", "C18_neg")])

hmdb <- function (x) {
  x <- x[-c(1:4),-c(1:4,6:7)]
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- x[complete.cases(x[,1]), ]
  x[,-1] <- sapply(x[,-1],as.numeric)
  x <- as.data.frame(t(x))
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- tibble::rownames_to_column(x, "samples")
  x <- x %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
  x <- merge(x, outcome, by = "samples", all.x = TRUE)
  x <- x %>% relocate(status, .after = samples)
  x <- x[complete.cases(x[,2]), ]
  x <- x %>% 
    mutate(status = replace(status, status == 1, "Converter"))
  x <- x %>% 
    mutate(status = replace(status, status == 2, "Persistently_uninfected"))
}





pool <- function(x) {
  x <- x[-c(1:4),-c(1:4,6:7)]
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- x[complete.cases(x[,1]), ]
  x <- as.data.frame(t(x))
  x <-`colnames<-`(x,x[1,])
  x <- x[-c(1),]
  x <- tibble::rownames_to_column(x, "samples")
  x <- x %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
  x <- x %>% filter(str_detect(samples,"^PREFB"))
} 

HMDB_pool <- parSapply(cl,HMDB_ID,pool)
x <- HMDB_ID$HIL_neg[-c(1:4),-c(1:4,6:7)]
x <-`colnames<-`(x,x[1,])
x <- x[-c(1),]
x <- x[complete.cases(x[,1]), ]
x$HMDB_ID <- gsub("\\*","",x$HMDB_ID)
x <- as.data.frame(t(x))
x <-`colnames<-`(x,x[1,])
x <- x[-c(1),]
x <- tibble::rownames_to_column(x, "samples")
x <- x %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
HMDB_pool$HIL_neg <- x %>% filter(str_detect(samples,"^PREFA"))

HMDB_ID <- parSapply(cl,HMDB_ID,hmdb)
colnames(HMDB_ID$HIL_pos) <- gsub("\\*","",colnames(HMDB_ID$HIL_pos))
colnames(HMDB_ID$HIL_neg) <- gsub("\\*","",colnames(HMDB_ID$HIL_neg))
colnames(HMDB_ID$C8_pos) <- gsub("\\*","",colnames(HMDB_ID$C8_pos))
colnames(HMDB_ID$C18_neg) <- gsub("\\*","",colnames(HMDB_ID$C18_neg))

#HMDB_ID$C8_pos <- HMDB_ID$C8_pos %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
#HMDB_ID$HIL_pos <- HMDB_ID$HIL_pos %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
#HMDB_ID$HIL_neg <- HMDB_ID$HIL_neg %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))
#HMDB_ID$C18_neg <- HMDB_ID$C18_neg %>% select(!starts_with(c("redundant","Internal","NA","REMOVE")))

HMDB_ID$HIL_pos$status  <- as.factor(HMDB_ID$HIL_pos$status)
HMDB_ID$HIL_neg$status  <- as.factor(HMDB_ID$HIL_neg$status)
HMDB_ID$C8_pos$status   <- as.factor(HMDB_ID$C8_pos$status)
HMDB_ID$C18_neg$status  <- as.factor(HMDB_ID$C18_neg$status)
HMDB_ID$HIL_pos[,-c(1:2)] <- parSapply(cl,HMDB_ID$HIL_pos[,-c(1:2)],as.numeric)
HMDB_ID$HIL_neg[,-c(1:2)] <- parSapply(cl,HMDB_ID$HIL_neg[,-c(1:2)],as.numeric)
HMDB_ID$C8_pos[,-c(1:2)]  <- parSapply(cl,HMDB_ID$C8_pos[,-c(1:2)],as.numeric)
HMDB_ID$C18_neg[,-c(1:2)] <- parSapply(cl,HMDB_ID$C18_neg[,-c(1:2)],as.numeric)


# # 4 finding doubles----------------------


HIL_pos$method <- "HILpos"
HIL_neg$method <- "HILneg"
C8_pos$method <- "C8pos"
C18_neg$method <- "C18neg"

all_hmdb <- add_row(HIL_pos,HIL_neg)
all_hmdb <- add_row(all_hmdb,C18_neg)
all_hmdb <- add_row(all_hmdb,C8_pos)

duplicate <- all_hmdb$HMDB_ID[duplicated(all_hmdb$HMDB_ID)]

x <- all_hmdb %>% filter(all_hmdb$HMDB_ID %in% duplicate)
x.hilpos <- x %>% filter(x$method == "HILpos")
x.hilneg <- x %>% filter(x$method == "HILneg")
x.c18neg <- x %>% filter(x$method == "C18neg")
x.c8pos <- x %>% filter(x$method == "C8pos")

HIL_pos_pool <- HMDB_pool$HIL_pos
HIL_neg_pool <- HMDB_pool$HIL_neg
C18_neg_pool <- HMDB_pool$C18_neg
C8_pos_pool <- HMDB_pool$C8_pos

HIL_pos_pool <- HIL_pos_pool %>% select(c(x.hilpos$HMDB_ID))
HIL_pos_pool <- as.data.frame(sapply(HIL_pos_pool,as.numeric))
HIL_neg_pool <- HIL_neg_pool %>% select(c(x.hilneg$HMDB_ID))
HIL_neg_pool <- as.data.frame(sapply(HIL_neg_pool,as.numeric))
C18_neg_pool <- C18_neg_pool %>% select(c(x.c18neg$HMDB_ID))
C18_neg_pool <- as.data.frame(sapply(C18_neg_pool,as.numeric))
C8_pos_pool <- C8_pos_pool %>% select(c(x.c8pos$HMDB_ID))
C8_pos_pool <- as.data.frame(sapply(C8_pos_pool,as.numeric))


HIL_pos_pool <- as.data.frame(apply(HIL_pos_pool, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100}))
HIL_pos_pool$method <- "HILpos"
HIL_neg_pool <- as.data.frame(apply(HIL_neg_pool, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100}))
HIL_neg_pool$method <- "HILneg"
C8_pos_pool <- as.data.frame(apply(C8_pos_pool, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100}))
C8_pos_pool$method <- "C8pos"
C18_neg_pool <- as.data.frame(apply(C18_neg_pool, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)*100}))
C18_neg_pool$method <- "C18neg"


names(HIL_pos_pool)[names(HIL_pos_pool) == "apply(HIL_pos_pool, 2, function(x) {     sd(x, na.rm = T)/mean(x, na.rm = T) * 100 })"] <- "cv"
names(HIL_neg_pool)[names(HIL_neg_pool) == "apply(HIL_neg_pool, 2, function(x) {     sd(x, na.rm = T)/mean(x, na.rm = T) * 100 })"] <- "cv"
names(C8_pos_pool)[names(C8_pos_pool) == "apply(C8_pos_pool, 2, function(x) {     sd(x, na.rm = T)/mean(x, na.rm = T) * 100 })"] <- "cv"
names(C18_neg_pool)[names(C18_neg_pool) == "apply(C18_neg_pool, 2, function(x) {     sd(x, na.rm = T)/mean(x, na.rm = T) * 100 })"] <- "cv"

HIL_pos_pool <- tibble::rownames_to_column(HIL_pos_pool, "sample")
HIL_neg_pool <- tibble::rownames_to_column(HIL_neg_pool, "sample")
C18_neg_pool <- tibble::rownames_to_column(C18_neg_pool, "sample")
C8_pos_pool <- tibble::rownames_to_column(C8_pos_pool, "sample")



all_hmdb_pool <- add_row(HIL_pos_pool,HIL_neg_pool)
all_hmdb_pool <- add_row(all_hmdb_pool,C18_neg_pool)
all_hmdb_pool <- add_row(all_hmdb_pool,C8_pos_pool)
all_hmdb_pool<- all_hmdb_pool%>% arrange(desc(all_hmdb_pool$sample))
str(all_hmdb_pool)
all_hmdb_pool$sample <- as.factor(all_hmdb_pool$sample)

selected <- all_hmdb_pool %>% 
  group_by(sample) %>% 
  slice(which.min(cv))
selected$action <- 1

all_hmdb_pool <- left_join(all_hmdb_pool,selected)
all_hmdb_pool$action[is.na(all_hmdb_pool$action)] = 0

### all_hmdb_pool modified outside in excel###
write.csv (all_hmdb_pool,paste0(dir,"duplicate.csv"))
### all_hmdb_pool modified outside in excel###

HILpos_discard <- all_hmdb_pool %>% filter(all_hmdb_pool$action==0 & all_hmdb_pool$method == "HILpos")
HILneg_discard <- all_hmdb_pool %>% filter(all_hmdb_pool$action==0 & all_hmdb_pool$method == "HILneg")
C8pos_discard <- all_hmdb_pool %>% filter(all_hmdb_pool$action==0 & all_hmdb_pool$method == "C8pos")
C18neg_discard <- all_hmdb_pool %>% filter(all_hmdb_pool$action==0 & all_hmdb_pool$method == "C18neg")

# # 3.1. save metabolites names for annotated and discard duplicates -------------------------------------------------

HMDB_ID$HIL_pos <- HMDB_ID$HIL_pos %>% select(!c(HILpos_discard$sample))
HMDB_ID$HIL_neg <- HMDB_ID$HIL_neg %>% select(!c(HILneg_discard$sample))
HMDB_ID$C8_pos  <- HMDB_ID$C8_pos  %>% select(!c(C8pos_discard$sample))
HMDB_ID$C18_neg <- HMDB_ID$C18_neg %>% select(!c(C18neg_discard$sample))


write.csv(HMDB_ID$HIL_pos,paste0(dir,"1. HILIC-pos_annotated.HMDB.csv"), row.names = FALSE)
write.csv(HMDB_ID$HIL_neg,paste0(dir,"1. HILIC-neg_annotated.HMDB.csv"), row.names = FALSE)
write.csv(HMDB_ID$C8_pos,paste0(dir,"1. C8-pos_annotated.HMDB.csv"), row.names = FALSE)
write.csv(HMDB_ID$C18_neg,paste0(dir,"1. C18-neg_annotated.HMDB.csv"), row.names = FALSE)
stopCluster(cl)

# QC drift histogram -----------------------------------------------------
sim_HM_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/2
    x
  })
  return(result)
}



#HILIC-positive
HIL_pos <- HIL_pos[-c(1:2,4:5),-c(1:6)]
HIL_pos <-`colnames<-`(HIL_pos,HIL_pos[1,])
HIL_pos <- HIL_pos[-c(1),]
HIL_pos <- HIL_pos[complete.cases(HIL_pos[,1]), ]
HIL_pos[,-1] <- sapply(HIL_pos[,-1],as.numeric)
str(HIL_pos)
HIL_pos$Sample_type <- as.factor(HIL_pos$Sample_type)
log_HIL_pos <- log(HIL_pos[-c(1)])
log_HIL_pos <- sim_HM_wrapper(log_HIL_pos)

myColors <- ifelse(colnames(log_HIL_pos)=="QC_drift_correction" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(log_HIL_pos)=="QC_pooled_plasma", rgb(0.8,0.1,0.3,0.6),
                          "grey90" ) )
# Build the plot
png(file=paste0(dir,"histogram_QC_HIL_pos.png"),width=1200, height=900)
boxplot(log_HIL_pos, 
        col=myColors , 
        ylab="Metabolites" , xlab="- Injection Order -")
dev.off()


#HILIC-negative
HIL_neg <- HIL_neg[-c(1:2,4:5),-c(1:6)]
HIL_neg <-`colnames<-`(HIL_neg,HIL_neg[1,])
HIL_neg <- HIL_neg[-c(1),]
HIL_neg <- HIL_neg[complete.cases(HIL_neg[,1]), ]
HIL_neg[,-1] <- sapply(HIL_neg[,-1],as.numeric)
str(HIL_neg)
HIL_neg$Sample_type <- as.factor(HIL_neg$Sample_type)
log_HIL_neg <- log(HIL_neg[-c(1)])
log_HIL_neg <- sim_HM_wrapper(log_HIL_neg)
myColors <- ifelse(colnames(log_HIL_neg)=="QC_drift_correction" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(log_HIL_neg)=="QC_pooled_plasma", rgb(0.8,0.1,0.3,0.6),
                          "grey90" ) )
# Build the plot
png(file=paste0(dir,"histogram_QC_HIL_neg.png"),width=1200, height=900)
boxplot(log_HIL_neg, 
        col=myColors , 
        ylab="Metabolites" , xlab="- Injection Order -")
dev.off()


###
#C8-positive
C8_pos <- C8_pos[-c(1:2,4:5),-c(1:6)]
C8_pos <-`colnames<-`(C8_pos,C8_pos[1,])
C8_pos <- C8_pos[-c(1),]
C8_pos <- C8_pos[complete.cases(C8_pos[,1]), ]
C8_pos[,-1] <- sapply(C8_pos[,-1],as.numeric)
str(C8_pos)
C8_pos$Sample_type <- as.factor(C8_pos$Sample_type)
log_C8_pos <- log(C8_pos[-c(1)])
log_C8_pos <- sim_HM_wrapper(log_C8_pos)
myColors <- ifelse(colnames(log_C8_pos)=="QC_drift_correction" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(log_C8_pos)=="QC_pooled_plasma", rgb(0.8,0.1,0.3,0.6),
                          "grey90" ) )
# Build the plot
png(file=paste0(dir,"histogram_QC_C8_pos.png"),width=1200, height=900)
boxplot(log_C8_pos, 
        col=myColors , 
        ylab="Metabolites" , xlab="- Injection Order -")
dev.off()




###
#C18-negative
C18_neg <- C18_neg[-c(1:2,4:5),-c(1:6)]
C18_neg <-`colnames<-`(C18_neg,C18_neg[1,])
C18_neg <- C18_neg[-c(1),]
C18_neg <- C18_neg[complete.cases(C18_neg[,1]), ]
C18_neg[,-1] <- sapply(C18_neg[,-1],as.numeric)
str(C18_neg)
C18_neg$Sample_type <- as.factor(C18_neg$Sample_type)
log_C18_neg <- log(C18_neg[-c(1)])
log_C18_neg <- sim_HM_wrapper(log_C18_neg)
myColors <- ifelse(colnames(log_C18_neg)=="QC_drift_correction" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(log_C18_neg)=="QC_pooled_plasma", rgb(0.8,0.1,0.3,0.6),
                          "grey90" ) )
# Build the plot
png(file=paste0(dir,"histogram_QC_C18_neg.png"),width=1200, height=900)
boxplot(log_C18_neg, 
        col=myColors , 
        ylab="Metabolites" , xlab="- Injection Order -")

dev.off()

#####################PCA#################

#HIL_pos
t_log_HIL_pos <- as.data.frame(t(log_HIL_pos))
t_log_HIL_pos$groups <- colnames(log_HIL_pos)
t_log_HIL_pos <- t_log_HIL_pos %>% relocate(groups,.before = V1)
str(t_log_HIL_pos)
pca_res <- prcomp(t_log_HIL_pos[-c(1)], scale. = TRUE)
png(file=paste0(dir,"PCA_QC_HIL_pos.png"),width=600, height=400)
autoplot(pca_res,data=t_log_HIL_pos,colour='groups')
dev.off()

#HIL_neg
t_log_HIL_neg <- as.data.frame(t(log_HIL_neg))
t_log_HIL_neg$groups <- colnames(log_HIL_neg)
t_log_HIL_neg <- t_log_HIL_neg %>% relocate(groups,.before = V1)
str(t_log_HIL_neg)
pca_res <- prcomp(t_log_HIL_neg[-c(1)], scale. = TRUE)
png(file=paste0(dir,"PCA_QC_HIL_neg.png"),width=600, height=400)
autoplot(pca_res,data=t_log_HIL_neg,colour='groups')
dev.off()

#C8_pos
t_log_C8_pos <- as.data.frame(t(log_C8_pos))
t_log_C8_pos$groups <- colnames(log_C8_pos)
t_log_C8_pos <- t_log_C8_pos %>% relocate(groups,.before = V1)
str(t_log_C8_pos)
pca_res <- prcomp(t_log_C8_pos[-c(1)], scale. = TRUE)
png(file=paste0(dir,"PCA_QC_C8_pos.png"),width=600, height=400)
autoplot(pca_res,data=t_log_C8_pos,colour='groups')
dev.off()

#C18_neg
C18_neg <- "1. 20_0401_INFECT_C18-neg_metabolomics_results.xlsx"
C18_neg <- read_excel(paste0(dir,C18_neg))
C18_neg <- C18_neg[-c(1:2,4),-c(1:6)]
C18_neg <-`colnames<-`(C18_neg,C18_neg[2,])
C18_neg <- C18_neg[-c(2),]
C18_neg <- C18_neg[complete.cases(C18_neg[,1]), ]
t_C18_neg <- as.data.frame(t(C18_neg))
t_C18_neg <-`colnames<-`(t_C18_neg,t_C18_neg[1,])
t_C18_neg <- t_C18_neg[-c(1),]
t_C18_neg[,-1] <- sapply(t_C18_neg[,-1],as.numeric)
str(t_C18_neg)
t_C18_neg[-c(1)] <- sim_HM_wrapper(t_C18_neg[-c(1)])
t_C18_neg[-c(1)] <- log(t_C18_neg[-c(1)])
png(file=paste0(dir,"PCA_QC_C18_neg.png"),width=600, height=400)
pca_res <- prcomp(t_C18_neg[-c(1)], scale. = TRUE)
autoplot(pca_res,data=t_C18_neg,colour='Sample_type')
dev.off()



t_C18_neg_noQC <- t_C18_neg[which(t_C18_neg$Sample_type =="Study_sample-INFECT"),]


#kmeans clustering
#C18
pacman::p_load(cluster,factoextra)
k2_2_clusters <- kmeans(t_C18_neg_noQC[-c(1)], centers = 2, nstart = 25)
k2_3_clusters <- kmeans(t_C18_neg_noQC[-c(1)], centers = 3, nstart = 25)
C18_k2_2cluster <- as.data.frame(k2_2_clusters$cluster)
C18_k2_2cluster <- tibble::rownames_to_column(C18_k2_2cluster, "IdInfect")
C18_k2_3cluster <- as.data.frame(k2_3_clusters$cluster)
C18_k2_3cluster <- tibble::rownames_to_column(C18_k2_3cluster, "IdInfect")
C18_k2_cluster <- merge(C18_k2_2cluster,C18_k2_3cluster,by="IdInfect")

png(file=paste0(dir,"k2_C18_neg.png"),width=600, height=400)
fviz_cluster(k2_2_clusters, data = t_C18_neg_noQC[-c(1)])
dev.off()

kmetabolite_2cluster <- as.data.frame(k2_2_clusters[["centers"]])
str(metabolite_2cluster)
write.csv(metabolite_2cluster,paste0(dir,"metabolite2cluster.csv"))
TableOne_cluster <- CreateTableOne(vars = colnames(metabolite_2cluster),data=metabolite_2cluster)
Stat_cluster <-print(TableOne_cluster, nonnormal = colnames(metabolite_2cluster)) %>% as.data.frame()


t_C18_neg_noQC <- tibble::rownames_to_column(t_C18_neg_noQC, "IdInfect")
t_C18_neg_noQC <- merge(t_C18_neg_noQC,C18_k2_2cluster,by="IdInfect")
t_C18_neg_noQC <- tibble::column_to_rownames(t_C18_neg_noQC,"IdInfect")
t_C18_neg_noQC <- as.data.frame(t(t_C18_neg_noQC))
t_C18_neg_noQC <-`colnames<-`(t_C18_neg_noQC,t_C18_neg_noQC[135,])
t_C18_neg_noQC <- t_C18_neg_noQC[-c(135),]
t_C18_neg_noQC <- t_C18_neg_noQC[-c(1),]

myColors2 <- ifelse(colnames(log_C18_neg)=="1" , rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(log_C18_neg)=="2", rgb(0.8,0.1,0.3,0.6),
                          "grey90"))

str(t_C18_neg_noQC)
t_C18_neg_noQC <- sapply(t_C18_neg_noQC,as.numeric)
# Build the plot
col_order <- c("1", "2")
my_data2 <- t_C18_neg_noQC[ ,order(colnames(t_C18_neg_noQC))]
my_data2

boxplot(my_data2,
        col=myColors2 , 
        ylab="Metabolites" , xlab="- Injection Order -")

write.csv(k2_2_clusters[["cluster"]],paste0(dir,"cluster_ID.csv"))
C18 <- read.csv(paste0(dir,"1. C18-neg_annotated.csv"))
cluster <- as.data.frame(k2_2_clusters[["cluster"]])
cluster <- tibble::rownames_to_column(cluster,"samples")

C18 <- merge(C18,cluster, by="samples")
colnames(C18)
p_load(naniar,ggplot2)
gg_miss_fct(x = C18[3:135], fct = C18$`k2_2_clusters[["cluster"]]`)
vis_miss(C18[order(C18$`k2_2_clusters[["cluster"]]`),])


x<- C18 %>%
  group_by(`k2_2_clusters[["cluster"]]`) %>%
  miss_var_summary()            


HILIC <- read.csv(paste0(dir,"1. C18-neg_raw.csv"))
?vis_miss
gg_miss_fct(x = HILIC, fct = C18$`k2_2_clusters[["cluster"]]`)
vis_miss(HILIC)
?vis_miss
x <- miss_var_summary(HILIC)
