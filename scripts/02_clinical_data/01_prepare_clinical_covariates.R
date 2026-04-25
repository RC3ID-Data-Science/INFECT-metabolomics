
# # load packages ---------------------------------------------------------
library(pacman)
pacman::p_load(tidyverse, readxl, psych, devtools, openxlsx, limma)
pacman::p_load(dlookr,dplyr, devtools,openxlsx,data.table,tableone,tidyverse,janitor)

# # File name and locations -----------------------------------------------
dir            <- "C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\"

all            <- "2. Eligible cases and contacts n1347 - repaired.xlsx" 
outcome        <- "2. assigned outcome.xlsx"
exposure_score <- "2. exposure score.xlsx"
contact.CRF    <- "2. F9 contact baseline CRF.xlsx"
index.CRF      <- "2. F3 CRF index case.xlsx" 
microbiology   <- "2. microbiology results bandung.xlsx" 
index.xray     <- "2. F4 index xray findings.xlsx"
vit.D          <- "2. Vitamin D.xlsx"
iron           <- "2. serum iron and inflam_2.xlsx"
micro          <- "2. micronutrients.xlsx"
wbsa           <- "2. wbsa wholeblood baseline.xlsx"
wgs            <- "2. wgs lineages index and exclusions.xlsx"
linkage        <- "2. linkage index contact.xlsx"
hematology     <- "2. haematology.xlsx"
socioenvironment <- "2. socioenvironment.xlsx"
IGRA <- "IGRA_INFECT_strict_cutoff.csv"

# (2) data klinis
list2env(sapply(list(all            = paste0(dir,all), 
                     outcome        = paste0(dir,outcome),
                     exposure_score = paste0(dir,exposure_score),
                     contact.CRF    = paste0(dir,contact.CRF),
                     index.CRF      = paste0(dir,index.CRF),
                     microbiology   = paste0(dir,microbiology),
                     index.xray     = paste0(dir,index.xray),
                     vit.D          = paste0(dir,vit.D),
                     iron           = paste0(dir,iron),
                     micro          = paste0(dir,micro),
                     wbsa           = paste0(dir,wbsa),
                     wgs            = paste0(dir,wgs),
                     linkage        = paste0(dir,linkage),
                     hematology     = paste0(dir,hematology),
                  socioenvironment  = paste0(dir,socioenvironment)),
                read_xlsx),envir=.GlobalEnv)

IGRA <- read.csv(paste0(dir,IGRA))
IGRA$deltaIGRA <- IGRA$follow_up-IGRA$baseline
HILIC <- read.csv(paste0(dir,"1. HILIC-pos_annotated.csv"))
metabolomic_case <- HILIC[,1:2]
metabolomic_case$status <- NULL

# # Adding variables --------------------------------------------------------
#adding ageindex, vitd, iron --- add more here!!
ageindex <- c("IdPatient","ageindex")
ageindex <- select(index.CRF,matches(ageindex))
index.xray$cavity <- ifelse(index.xray$f4_inf6n=="1", "1", ifelse(index.xray$f4_inf7n=="1","1",ifelse(index.xray$f4_inf8n=="1","1",ifelse(index.xray$f4_inf9n=="1","1",ifelse(index.xray$f4_inf10n=="1","1",ifelse(index.xray$f4_inf11=="1","1","0"))))))
index.xray <- index.xray[,c(1,29)]
index.xray$cavity<-as.numeric(index.xray$cavity)


wgs <- wgs[,c(1,3)]
all <- left_join(x=all,y=ageindex, by = c("linkage_index_contact.IdPatient" = "IdPatient"))
all <- left_join(x=all,y=vit.D, by = c("linkage_index_contact.IdInfect" = "IdInfect"))
all <- left_join(x=all,y=iron, by = c("linkage_index_contact.IdInfect" = "IdInfect"))
all <- left_join(x=all,y=wgs, by = c("linkage_index_contact.IdPatient" = "IdIndex"))
all <- left_join(x=all,y=index.xray, by = c("linkage_index_contact.IdPatient" = "IdPatient"))
all <- left_join(x=all,y=hematology, by = c("linkage_index_contact.IdInfect" = "IdInfect"))
all <- left_join(x=all,y=socioenvironment, by = c("linkage_index_contact.IdPatient" = "IdPatient"))
all <- left_join(x=all,y=wbsa, by = c("linkage_index_contact.IdInfect" = "IdInfect"))

#change variable names
colnames(all)
variable.metabolomic <- c("linkage_index_contact.IdInfect","outcome2","pred14","f9_inf_sex","agecontact","f9_inf17","f9_inf14","f9_inf9","f9_inf10","f9_inf29","f9_inf32","f9_inf34","f9_inf35","ageindex","inf_sex","Field0","cavity","f4_inf12","25 OHD3","25 OHD2","FER_ug/L","sTFR_as_Ramco_assay_equivalents_mg/L","Body_Iron_Stores_BIS_mg/kg_body_weight","RBP_as_retinol_equivalents_umol/L","CRP_mg/L","AGP_g/L", "Name","hemoglobin","hematocrit","leucocytes","platelet","neutrophils","lymphocytes","monocyte","f5_inf17","f5_inf18","f5_inf19","f5_inf20","f5_inf33","f5_inf34","f5_inf35","f5_inf31","ifnnil","ifnbcg","ifnmtb","ifnstrep","ifnecoli","ifncan","tnfnil","tnfbcg","tnfmtb","tnfstrep","tnfecoli","tnfcan","il1bnil","il1bbcg","il1bmtb","il1bstrep","il1becoli","il1bcan","il1ranil","il1rabcg","il1ramtb","il1rastrep","il1raecoli","il1racan","il1anil","il1abcg","il1amtb","il1astrep","il1aecoli","il1acan","il6nil","il6bcg","il6mtb","il6strep","il6ecoli","il6can","il8nil","il8bcg","il8mtb","il8strep","il8ecoli","il8can","il10nil","il10bcg","il10mtb","il10strep","il10ecoli","il10can","il32nil","il32bcg","il32mtb","il32pneumo","il32ecoli","il32candida","nobbcg","nobmtb","nobstrep","nobecoli","linkage_index_contact.IdPatient")
variable.metabolomic2 <- c("IdInfect","status","exp_score","sex_contact","agecontact","BCG_scar","smoking","sleep_prox","hours_with_caseindex","weight","height","random_blood_glucose","HbA1c","ageindex","sex_index","smear_result_3","cavity","extent","25.OHD3","25.OHD2","FER_ug/L","sTFR_as_Ramco_assay_equivalents_mg/L","Body_Iron_Stores_BIS_mg/kg_body_weight","RBP_as_retinol_equivalents_umol/L","CRP_mg/L","AGP_g/L", "Mtb strain","hemoglobin","hematocrit","leucocytes","platelet","neutrophils","lymphocytes","monocyte","house_area","no.bedrooms","no.rooms","windows","cooking","lightsource","watersource","household_income","ifnnil","ifnbcg","ifnmtb","ifnstrep","ifnecoli","ifncan","tnfnil","tnfbcg","tnfmtb","tnfstrep","tnfecoli","tnfcan","il1bnil","il1bbcg","il1bmtb","il1bstrep","il1becoli","il1bcan","il1ranil","il1rabcg","il1ramtb","il1rastrep","il1raecoli","il1racan","il1anil","il1abcg","il1amtb","il1astrep","il1aecoli","il1acan","il6nil","il6bcg","il6mtb","il6strep","il6ecoli","il6can","il8nil","il8bcg","il8mtb","il8strep","il8ecoli","il8can","il10nil","il10bcg","il10mtb","il10strep","il10ecoli","il10can","il32nil","il32bcg","il32mtb","il32pneumo","il32ecoli","il32candida","nobbcg","nobmtb","nobstrep","nobecoli","IdPatient")
metabolomic <- select(all,matches(variable.metabolomic, ignore.case = FALSE))
metabolomic <- metabolomic %>%
  rename_at(vars(variable.metabolomic), ~(variable.metabolomic2))

all <- select(all,matches(variable.metabolomic, ignore.case = FALSE))
all <- all %>%
  rename_at(vars(variable.metabolomic), ~(variable.metabolomic2))


#adding BMI and Diabetes (variables that need to be calculated or with IF function add here!!)
metabolomic <- mutate(metabolomic,BMI_contact=metabolomic$weight/(metabolomic$height^2))
metabolomic<-mutate(metabolomic,DM=if_else(metabolomic$HbA1c>6.4,2,if_else(metabolomic$HbA1c<5.7,0,1),missing=if_else(metabolomic$random_blood_glucose<101,0,NA)))
#2=diabetes,1=prediabetes,0=no diabetes#
all <- mutate(all,BMI_contact=all$weight/(all$height^2))
all <-mutate(all,DM=if_else(all$HbA1c>6.4,2,if_else(all$HbA1c<5.7,0,1),missing=if_else(all$random_blood_glucose<101,0,NA)))

all$agegrp <-cut(all$agecontact, c(4,8,17,29,39,49,101),
        labels=c("5-9", "9-18", "18-30","30-40","40-50","50+"))


# Choose cases with converter and persistenly uninfected outcome --------


#all$IdPatient <- NULL
all$random_blood_glucose <- NULL
all$HbA1c <- NULL
all$height <- NULL
all$weight <- NULL
all$f5_inf33b <- NULL
all$f5_inf34b <- NULL
all$f5_inf35b <- NULL
all$lightsource <- NULL

# # Change value names ----------------------------------------------------
test<-all$status =="1"
all$status[test] <- "Converter"
test<-all$status =="2"
all$status[test] <- "Persistently_uninfected"
test<-all$status =="0"
all$status[test] <- "Baseline_positive"

test<-all$sex_contact =="0"
all$sex_contact[test] <- "Female"
test<-all$sex_contact =="1"
all$sex_contact[test] <- "Male"

test<-all$sex_index =="0"
all$sex_index[test] <- "Female"
test<-all$sex_index =="1"
all$sex_index[test] <- "Male"


test<-all$BCG_scar =="2"
all$BCG_scar[test] <- NA
test<-all$BCG_scar =="1"
all$BCG_scar[test] <- "2"
test<-all$BCG_scar =="0"
all$BCG_scar[test] <- "1"



test<-all$smoking =="1"
all$smoking[test] <- "3"
test<-all$smoking =="2"
all$smoking[test] <- "2"
test<-all$smoking =="0"
all$smoking[test] <- "1"


test<-all$sleep_prox =="1"
all$sleep_prox[test] <- "2"
test<-all$sleep_prox =="2"
all$sleep_prox[test] <- "2"
test<-all$sleep_prox =="0"
all$sleep_prox[test] <- "1"

test<-all$cavity %in% c("0","2","3","4","5")
all$cavity[test] <- "0"
test<-all$cavity =="1"
all$cavity[test] <- "2"
test<-all$cavity == "0"
all$cavity[test] <- "1"

test<-all$DM =="2"
all$DM[test] <- "3"
test<-all$DM =="1"
all$DM[test] <- "2"
test<-all$DM =="0"
all$DM[test] <- "1"


test <- all$smear_result_3 == "Positive_3"
all$smear_result_3[test] <- "5"
test <- all$smear_result_3 == "Positive_2"
all$smear_result_3[test] <- "4"
test <- all$smear_result_3 == "Positive_1"
all$smear_result_3[test] <- "3"
test <- all$smear_result_3 == "Scanty/Positive_CSF"
all$smear_result_3[test] <- "2"
test <- all$smear_result_3 == "Negative"
all$smear_result_3[test] <- "1"

test <- all$cooking == "1"
all$cooking[test] <- "Gas"
test <- all$cooking == "4"
all$cooking[test] <- "Wood"

test <- all$watersource == "1"
all$watersource[test] <- "PDAM"
test <- all$watersource == "2"
all$watersource[test] <- "Bottled_water"
test <- all$watersource == "3"
all$watersource[test] <- "well"
test <- all$watersource == "5"
all$watersource[test] <- "Other"

all$status <-
  if_else(all$status== 4,"Indeterminate",
          if_else(all$status == 5,"Active TB",
                  if_else(all$status == 6,"Unevaluated symptoms",
                          if_else(all$status == 9,"Lost to follow-up",
                                  all$status)
                  )
          )
          
  )


colnames(all)
write.csv(all,paste0(dir,"all.csv"))
persistently_uninfected <- filter(all,status=="Persistently_uninfected")
converter <- filter(all,status=="Converter")
baseline_positive <- all %>% filter(status=="Baseline_positive")
write.csv(baseline_positive,paste0(dir,"klinis_baseline_positive.csv"))
all_2 <- all
all<-merge(converter,persistently_uninfected,all = TRUE)

# # change data types -----------------------------------------------------
str(all)
colnames(all)
all$cavity <- as.factor(all$cavity)
all$IdInfect <- sapply(all$IdInfect,as.character)
all[,c(2,4,6,8,11,13,23,32,33)] <- lapply(all[,c(2,4,6,8,11,13,23,32,33)] ,as.factor)
all[,c(7,12,34,94,95)] <- lapply(all[,c(7,12,34,94,95)] ,ordered)



################# risk and exposure ################
colnames(all)
variable.metabolomic3 <- c("sex_contact","status","sex_index","agecontact","agegrp","exp_score","BCG_scar","smoking","BMI_contact","DM",
                           "sleep_prox","hours_with_caseindex","smear_result_3","cavity","extent","Mtb strain","platelet","neutrophils","lymphocytes","monocyte","agegrp")
#normality(all)
variable.metabolomic4 <- c(colnames(all[,-1]))
TableOne <- CreateTableOne(vars = variable.metabolomic4,strata="status",data=all)
Stat_metabolomic<-print(TableOne) %>% as.data.frame()

adults <- all %>% filter(all$agecontact > 18)
child <- all %>% filter(all$agecontact< 9 )
adults_male <- all %>% filter(all$agecontact > 18 & all$sex_contact=="Male")
adults_female <- all %>% filter(all$agecontact > 18 & all$sex_contact=="Female")

metabolomic_only <- all %>% filter(IdInfect %in% metabolomic_case$samples)

exposure_variable <- c("IdInfect","status","BMI_contact","agegrp","hours_with_caseindex","sex_index","household_income","cooking","no.bedrooms","no.rooms","windows","smear_result_3","cavity","extent","sleep_prox")

exposure_table <- select(all,matches(exposure_variable, ignore.case = FALSE))
exposure_table$total_rooms <- exposure_table$no.bedrooms + exposure_table$no.rooms
str(exposure_table)


write.csv(exposure_table,paste0(dir,"2. exposure_PCA.csv"), row.names = FALSE)
##############scoring risk ######################
exposure_table$score <- exposure_table$BMI_contact
#sex_index
exposure_table$score <-
  if_else(
    exposure_table$sex_index=="Male", exposure_table$score - 4, exposure_table$score - 0
  )
#agegrp
exposure_table$score <-
  if_else(
    exposure_table$agegrp=="9-18", exposure_table$score - 2, exposure_table$score - 0
  )
exposure_table$score <-
  if_else(
    exposure_table$agegrp=="50+", exposure_table$score - 2, exposure_table$score - 0
  )
#hours_with_case
exposure_table$score <-
  if_else(
    exposure_table$sleep_prox==2 & exposure_table$hours_with_caseindex>5, exposure_table$score - 4, exposure_table$score - 0
  )
#household_income
exposure_table$score <-
  if_else(
    exposure_table$household_income<3 , exposure_table$score - 2, exposure_table$score - 0
  )
#cooking
exposure_table$score <-
  if_else(
    exposure_table$cooking=="wood" , exposure_table$score - 2, exposure_table$score - 0
  )
#windows
exposure_table$score <-
  if_else(exposure_table$windows == 0, exposure_table$score - 3,
          if_else(exposure_table$total_rooms==2|3 & exposure_table$windows <=2,exposure_table$score - 3,
                  if_else(exposure_table$total_rooms==4 & exposure_table$windows <=3,exposure_table$score - 3,
                          if_else(exposure_table$total_rooms==5|6 & exposure_table$windows <=4,exposure_table$score - 3,exposure_table$score - 0))))

#risk_group
exposure_table$risk_group <- 
  if_else(
    exposure_table$score < 12, "high_risk",
    if_else(exposure_table$score >18, "low_risk", "medium_risk")
  )
############# exposure todia ###################

exposure_table$score_todi <- 0
#windows
exposure_table$score_todi <-
  if_else(exposure_table$windows == 0, exposure_table$score_todi +1,
          if_else(exposure_table$total_rooms==2|3 & exposure_table$windows <=2,exposure_table$score_todi + 1,
                  if_else(exposure_table$total_rooms==4 & exposure_table$windows <=3,exposure_table$score_todi + 1,
                          if_else(exposure_table$total_rooms==5|6 & exposure_table$windows <=4,exposure_table$score_todi +1 ,exposure_table$score_todi - 0))))
#hours with case
exposure_table$score_todi <-
  if_else(
    exposure_table$sleep_prox==2 & exposure_table$hours_with_caseindex>5, exposure_table$score_todi + 1, exposure_table$score_todi - 0
  )

#index infectivity
exposure_table$score_todi <-
  if_else(
    exposure_table$smear_result_3>=4 & exposure_table$cavity==2, exposure_table$score_todi + 1, exposure_table$score_todi - 0
  )


#exposure_group
exposure_table$exposure_group <- 
  if_else(
    exposure_table$score_todi >=  2, "high_exposure","low_exposure")

exposure_table_onlyrisk <- exposure_table[c("IdInfect","risk_group","exposure_group")]
all <- left_join(all,exposure_table_onlyrisk, by = c("IdInfect" = "IdInfect"))

variable.metabolomic5 <- c(colnames(exposure_table[,-1]))
variable.metabolomic4 <- c(colnames(all[,-1]))
#TableOne <- CreateTableOne(vars = colnames(),strata=c("status"),data=klinis)
#Stat_metabolomic<-print(TableOne) %>% as.data.frame()

#write.csv(Stat_metabolomic,paste0(dir,"table_all.csv"))

##########Mandalaka's score ################

#add IGRA results
IGRA$X <- NULL
IGRA_strict <- IGRA %>% select(IdInfect,strict_0.15_conversion,deltaIGRA,Log2FC_TBonly)
IGRA_strict$IdInfect <- as.character(IGRA_strict$IdInfect)
IGRA_strict <- arrange(IGRA_strict,IdInfect)
all <- left_join(all,IGRA_strict,by="IdInfect")

# # save clinical datasets ------------------------------------------------

sort <- c("IdInfect","status","exp_score","sex_contact","agecontact","BCG_scar","smoking","sleep_prox","hours_with_caseindex","ageindex","sex_index","smear_result_3","cavity","extent","25.OHD3","25.OHD2","FER_ug/L","sTFR_as_Ramco_assay_equivalents_mg/L","Body_Iron_Stores_BIS_mg/kg_body_weight","RBP_as_retinol_equivalents_umol/L","CRP_mg/L","AGP_g/L", "Mtb strain","hemoglobin","hematocrit","leucocytes","platelet","neutrophils","lymphocytes","monocyte","house_area","no.bedrooms","no.rooms","windows","cooking","watersource","household_income","ifnnil","ifnbcg","ifnmtb","ifnstrep","ifnecoli","ifncan","tnfnil","tnfbcg","tnfmtb","tnfstrep","tnfecoli","tnfcan","il1bnil","il1bbcg","il1bmtb","il1bstrep","il1becoli","il1bcan","il1ranil","il1rabcg","il1ramtb","il1rastrep","il1raecoli","il1racan","il1anil","il1abcg","il1amtb","il1astrep","il1aecoli","il1acan","il6nil","il6bcg","il6mtb","il6strep","il6ecoli","il6can","il8nil","il8bcg","il8mtb","il8strep","il8ecoli","il8can","il10nil","il10bcg","il10mtb","il10strep","il10ecoli","il10can","il32nil","il32bcg","il32mtb","il32pneumo","il32ecoli","il32candida","nobbcg","nobmtb","nobstrep","nobecoli","BMI_contact","DM","agegrp","risk_group","exposure_group","strict_0.15_conversion","deltaIGRA","Log2FC_TBonly","IdPatient")
all_sorted <- all %>% select(sort)
write.csv(all_sorted,paste0(dir,"2. all_klinis_subjekINFECT.csv"), row.names = FALSE)
write_rds(all_sorted,paste0(dir,"2. all_klinis_subjekINFECT.rds"))

write.csv(all_2,paste0(dir,"2. MEGA_all_klinis_subjekINFECT.csv"), row.names = FALSE)
write_rds(all_2,paste0(dir,"2. MEGA_all_klinis_subjekINFECT.rds"))

# # Descriptive analysis --------------------------------------------------

library(FactoMineR)
library(factoextra)

exposure_table$no.bedrooms <- NULL
exposure_table$no.rooms <- NULL
rownames(exposure_table) <- exposure_table$IdInfect
exposure_table$IdInfect <- NULL
exposure_table$status <- NULL
exposure_table$score <- NULL
exposure_table$risk_group <- NULL

str(exposure_table)
exposure_table <- lapply(exposure_table, as.numeric)
res.mca <- prcomp(na.omit(exposure_table), center = TRUE, scale = TRUE)

exposure_table$agegrp <- as.numeric(exposure_table$agegrp)
exposure_table$sex_index <- as.numeric(exposure_table$sex_index)
exposure_table$household_income <- as.numeric(exposure_table$household_income)
exposure_table$cooking <- as.numeric(exposure_table$cooking)
exposure_table$smear_result_3 <- as.numeric(exposure_table$smear_result_3)
exposure_table$cavity <- as.numeric(exposure_table$cavity)
exposure_table$sleep_prox <- as.numeric(exposure_table$sleep_prox)

fviz_eig(res.mca)



fviz_pca_ind(res.mca,axes = c(1,2),
             label = "none", # hide individual labels
)
fviz_pca_var(res.mca, axes = c(1,2))

# # K-mean clustering C18 negative-----------------------
metabolomic_C18_cluster <- merge(metabolomic,C18_k2_cluster,by="IdInfect")
colnames(metabolomic_C18_cluster)

#statistic k2 two clusters
variable.metabolomic4 <- c("status","sex_contact","agecontact","exp_score","BCG_scar","smoking","BMI_contact","DM",
                           "sleep_prox","hours_with_caseindex","smear_result_3","cavity","extent","Mtb strain")
TableOne_k2_2 <- CreateTableOne(vars = variable.metabolomic4,strata="k2_2_clusters$cluster",data=metabolomic_C18_cluster)
Stat_metabolomic_k2_2 <-print(TableOne_k2_2, nonnormal = variable.metabolomic4) %>% as.data.frame()

TableOne_k2_3 <- CreateTableOne(vars = variable.metabolomic4,strata="k2_3_clusters$cluster",data=metabolomic_C18_cluster)
Stat_metabolomic_k2_3 <-print(TableOne_k2_3, nonnormal = variable.metabolomic4) %>% as.data.frame()

