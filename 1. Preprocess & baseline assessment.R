############################################################
#Data pre-processing by study visit
#1.For individual metabolite analysis: 
#  - Inverse normal transformation by batch
#2.For creation of modules:
#  - FirstQuantil normalization across batches
#  - Inverse normal transofmration across batches
############################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# require('BiocManger')
#BiocManger::install("preprocessCore")
rm(list = ls())
library(preprocessCore)
setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
lipid <- read.csv("Data/Lipidomics_use.csv", header=T, sep=",") #1116*447
lipid <- lipid[!(lipid$Visit==1 & lipid$ID == '20477'),] #remove 1 outlier

####################### Func: Inverse Normal transform ##########################
invNorm = function(x){
  names(x) <- seq(1,length(x))  # Indexing each obs.
  x.new <- x[!is.na(x)]         # Separate missing from non-missing obs.
  x.na  <- x[is.na(x)]
  xnew.InvN <- qnorm(rank(x.new)/(length(x.new)+1))    # Perform INT on non-missing values
  x.InvN <- c(xnew.InvN,x.na)   # Combine INT-obs. with missing
  x.InvN <- x.InvN[order(as.numeric(names(x.InvN)))]  #sort by orginal order of the vector
  return(x.InvN)
}
####################### ----- LOD ----- ###########################
sum(apply(lipid[,6:425],2,class) != "integer") # 0 none of them are coded as non-interger
sum(apply(lipid[,6:425],2, function(x) sum(x < 0)),na.rm=T) #0 none of them are coded as negative value

##################### Multiple Imputation #########################
library(dplyr)
library(mice)
library(missForest)
library(foreign)

anno <- read.csv("../Data/Lipidomics_annotation_file.csv",header=T, sep=",")
anno$met_ID <- chartr(old = ".", new = "_", anno$lipid_identifier)

p_missing <- unlist(apply(lipid[,334:425],2,function(x) sum(is.na(x))/length(x)))
p_missing <- sort(p_missing[p_missing > 0], decreasing = T)

p_missing <- data.frame(met_ID = substring(names(p_missing),4),rate = p_missing)
p_missing <- merge(p_missing, anno[,c("met_ID","Annotation")],by.x = "met_ID")
p_missing <- p_missing[order(p_missing$rate,decreasing = T),]

missing_pattern <- md.pattern(lipid[,6:425])
missing_pattern[,ncol(missing_pattern)]

#### Impute with time and batch number
lipid$batch <- factor(lipid$batch)
imp <- missForest(lipid[,c(3,6:425,426)])
imp$OOBerror

#### replace the lipid measurements with imputed values
lipid.imp <- lipid
lipid.imp[,c(3,6:425)] <- imp$ximp[,1:421]
write.csv(lipid.imp,"Data/Lipidomics_use_imputed.csv",row.names = F) 

####################### Data processess ###########################
lipid <- read.csv("Data/Lipidomics_use_imputed.csv")
for (visit in c(0, 1, 2, 4)) {    # By visit
  # N = 318, 310, 209, 191
  DF <- lipid[lipid$Visit==visit,]
  
  ###1.Inverse-normal transformation (InvN)) by batch-----------------------------
  x.invN <- matrix(nrow=nrow(DF),ncol=length(6:425))
  colnames(x.invN) <- colnames(DF)[6:425]
  x.std <- matrix(nrow=nrow(DF),ncol=length(6:425))
  colnames(x.std) <- colnames(DF)[6:425]
  
  for (k in 1:3){                 # By batch (k)
    x.byBatch <- subset(DF,batch==k)[6:425]
    x.invN[DF$batch==k,] <- apply(x.byBatch,2,invNorm)
    x.std[DF$batch==k,] <- apply(x.byBatch,2,scale)
  }
  
  DF_InvN <- data.frame(ID = DF$ID, x.invN, wt = DF$wt, batch = DF$batch, GA_blood = DF$GA_blood)
  DF_Std <- data.frame(ID = DF$ID, x.std, wt = DF$wt, batch = DF$batch, GA_blood = DF$GA_blood)
  
  write.csv(DF_InvN,paste0("tmp/lipid_InvN",visit,".csv"), row.names=F) 
  write.csv(DF_Std,paste0("tmp/lipid_Std",visit,".csv"), row.names=F)
  
  ###2.Quantile normalization followed by InvN------------------------------------
  #Index <- c(6:333) #328 knonw lp_met
  Index <- c(6:425) #All 420 lipid metabolites
  #Index <- c(334:425) #All unknown lipid metabolites
  
  # Quantile normalization
  DF_T <- as.matrix(t(DF[,Index]))
  DF_QN <- t(normalize.quantiles(DF_T)) #missing values are left out automatically
  colnames(DF_QN) <- colnames(DF[,Index])
  
  # InvN
  DF_QN_InvN <- apply(DF_QN,2,invNorm)
  DF_QN_Std <- apply(DF_QN,2,scale)
  DF_QN_InvN <- data.frame(ID = DF$ID, DF_QN_InvN, wt = DF$wt, batch = DF$batch, GA_blood = DF$GA_blood)
  DF_QN_Std <- data.frame(ID = DF$ID, DF_QN_Std, wt = DF$wt, batch = DF$batch, GA_blood = DF$GA_blood)
  
  table(apply(DF_QN_InvN,2,function(vec){sum(is.na(vec))}))
  write.csv(DF_QN_InvN,paste0("tmp/lipid_QN_InvN",visit,".csv"), row.names=F)
  write.csv(DF_QN_Std, paste0("tmp/lipid_QN_Std", visit,".csv"), row.names = F)
}

###################### Baseline Characteristics #########################
library(haven)
library(coin)
library(spatstat)
library(weights)
library(stringr)
library(dplyr)
lipid <- read.csv("Data/Lipidomics_use_imputed.csv")
anthro <- read_sas("Data/formsingle.sas7bdat", col_select = c("SUBJECT_ID","ENRRACE_fm002","PreBMI","Age","Education",
                    "Paritycat","SMOKE6MO_fm019","Weightv0","Weightv1","Weightv2","Weightv3",
                    "Weightv4","Weightv5","Weightv6","NDSBWGT_fm024","Birthweightcat","Avg_NLENGTH_fm016",
                    "SITE_ID_fm002","Avg_NABDMC_fm016","Avg_NCHSTC_fm016","Avg_NHEADC_fm016","Avg_NMUACIR_fm016",
                    "Avg_NUMBCIR_fm016","SGA_Duryea","LGA_Duryea","INFSEXR1_fm024","GAv0","GAv1","GAv2","GAv4","ALCOL3M_fm011",
                    "DIET_DT_fm021_V0", "DIET_DT_fm021_V1", "DIET_DT_fm021_V2", "DIET_DT_fm021_V4",
                    "DIET_HR_fm021_V0", "DIET_HR_fm021_V1", "DIET_HR_fm021_V2", "DIET_HR_fm021_V4",
                    "DIET_MN_fm021_V0", "DIET_MN_fm021_V1", "DIET_MN_fm021_V2", "DIET_MN_fm021_V4",
                    "VIS_DT_fm021_V0", "VIS_DT_fm021_V1", "VIS_DT_fm021_V2", "VIS_DT_fm021_V4", 
                    "VISIT_HR_fm021_V0", "VISIT_HR_fm021_V1", "VISIT_HR_fm021_V2", "VISIT_HR_fm021_V4", 
                    "VISIT_MN_fm021_V0", "VISIT_MN_fm021_V1", "VISIT_MN_fm021_V2", "VISIT_MN_fm021_V4"))
names(anthro) <- c("SiteID", "ID", "Race","Age","Alcohol","Education","Weightv0","preBMI",
                   "Weightv1","Weightv2","Weightv3","Weightv4","Weightv5", "ab_circ",
                    "chest_circ", "head_circ", "Length", "arm_circ","um_circ", "Smoked", 
                   "vis0_dat", "vis1_dat", "vis2_dat", "vis4_dat", 
                   "vis0_hr", "vis1_hr", "vis2_hr", "vis4_hr",
                   "vis0_min","vis1_min","vis2_min","vis4_min",
                   "diet0_dat","diet1_dat","diet2_dat","diet4_dat",
                   "diet0_hr", "diet1_hr", "diet2_hr", "diet4_hr",
                   "diet0_min","diet1_min","diet2_min","diet4_min",
                   "bw", "Gender", "Weightv6", "bw_cat", "Parity", "GAv0","GAv1","GAv2","GAv4","SGA","LGA")
id_info <- read.csv("Data/pairID.csv")
diet0 <- as.POSIXlt(paste0(anthro$diet0_dat," ",str_pad(anthro$diet0_hr, 2, pad = "0"),":",str_pad(anthro$diet0_min,2,pad="0")), format = "%Y-%m-%d %H:%M")
visit0 <- as.POSIXlt(paste0(anthro$vis0_dat," ",str_pad(anthro$vis0_hr,2, pad="0"),":",str_pad(anthro$vis0_min, 2, pad = "0")), format = "%Y-%m-%d %H:%M")
sinceMealhr0 <- difftime(visit0, diet0,units="hours")
diet1 <- as.POSIXlt(paste0(anthro$diet1_dat," ",str_pad(anthro$diet1_hr, 2, pad = "0"),":",str_pad(anthro$diet1_min,2,pad="0")), format = "%Y-%m-%d %H:%M")
visit1 <- as.POSIXlt(paste0(anthro$vis1_dat," ",str_pad(anthro$vis1_hr,2, pad="0"),":",str_pad(anthro$vis1_min, 2, pad = "0")), format = "%Y-%m-%d %H:%M")
sinceMealhr1 <- difftime(visit1, diet1,units="hours")
diet2 <- as.POSIXlt(paste0(anthro$diet2_dat," ",str_pad(anthro$diet2_hr, 2, pad = "0"),":",str_pad(anthro$diet2_min,2,pad="0")), format = "%Y-%m-%d %H:%M")
visit2 <- as.POSIXlt(paste0(anthro$vis2_dat," ",str_pad(anthro$vis2_hr,2, pad="0"),":",str_pad(anthro$vis2_min, 2, pad = "0")), format = "%Y-%m-%d %H:%M")
sinceMealhr2 <- difftime(visit2, diet2,units="hours")
diet4 <- as.POSIXlt(paste0(anthro$diet4_dat," ",str_pad(anthro$diet4_hr, 2, pad = "0"),":",str_pad(anthro$diet4_min,2,pad="0")), format = "%Y-%m-%d %H:%M")
visit4 <- as.POSIXlt(paste0(anthro$vis4_dat," ",str_pad(anthro$vis4_hr,2, pad="0"),":",str_pad(anthro$vis4_min, 2, pad = "0")), format = "%Y-%m-%d %H:%M")
sinceMealhr4 <- difftime(visit4, diet4,units="hours")

anthro <- cbind(anthro,sinceMealhr0,sinceMealhr1,sinceMealhr2,sinceMealhr4)
anthro$Alcohol_bin <- factor(ifelse(anthro$Alcohol == 6,0, ifelse(anthro$Alcohol < 6,1,NA)))
anthro$Alcohol <- factor(anthro$Alcohol, ordered = T)
class(anthro$ID) <- "integer"
anthro <- anthro %>% 
  inner_join(id_info,by = "ID") 

anthro <- anthro %>% mutate(Nulliparious = Parity == 0,
                  LBW = bw_cat == 1, macro = bw_cat == 3)
anthro_add <- read_sas("Data/formsingle_anthro_20171115.sas7bdat", col_select = c("id","SMOKE6MO_fm019", "GWGv0","GWGv1","GWGv2","GWGv3","GWGv4","GWGv5","measuredate","SSF"))

names(anthro_add)[1:2] <- c("Smoked","ID")
anthro_add <- anthro_add %>%
  inner_join(id_info, by = "ID")

# order anthro and additional anthro by ID 
anthro <- anthro[order(anthro$ID),]
anthro_add <- anthro_add[order(anthro_add$ID),]
sum(anthro$ID == anthro_add$ID)

# merge the ordered anthro and anthro_add and only keep the useful covariates
anthro <- cbind(select(anthro, -c("Smoked","Weightv0","Weightv1","Weightv2","Weightv3","Weightv4","Weightv5","Weightv6")),
                anthro_add[,c("Smoked","GWGv0","GWGv1","GWGv2","GWGv4","measuredate","SSF")])
anthro$Smoked <- as.numeric(factor(anthro$Smoked == 1)) - 1
anthro$Education <- factor(anthro$Education, ordered = T)
anthro$Race <- factor(anthro$Race)


##### Table 1:-------------------------------------------------------------------#####
##### obtain the anthro variables for baseline assessment: age, nulliparious and smoked
#lipid <- read.csv("Data/Lipidomics_use_imputed.csv")

baseline_anthro <- lipid %>%
  filter(Visit == 0) %>%
  select(ID, age, wt, preBMI) %>%
  group_by(ID)

baseline_anthro <- baseline_anthro %>%
  inner_join(select(anthro, -c("Age","preBMI")),by="ID")

baseline_anthro$preBMIcat <- as.factor(baseline_anthro$preBMI < 30)
lBMI_idx <- baseline_anthro$preBMI < 30
hBMI_idx <- baseline_anthro$preBMI >= 30
sum(lBMI_idx)
sum(hBMI_idx)
## For age
  # median
  weighted.median(baseline_anthro$age,baseline_anthro$wt)
  weighted.median(baseline_anthro$age[lBMI_idx], baseline_anthro$wt[lBMI_idx])
  weighted.median(baseline_anthro$age[hBMI_idx], baseline_anthro$wt[hBMI_idx])
  # quantile
  weighted.quantile(baseline_anthro$age,baseline_anthro$wt, probs = c(0.25,0.75))
  weighted.quantile(baseline_anthro$age[lBMI_idx], baseline_anthro$wt[lBMI_idx], probs = c(0.25,0.75))
  weighted.quantile(baseline_anthro$age[hBMI_idx], baseline_anthro$wt[hBMI_idx], probs = c(0.25,0.75))
  # test
  wilcox_test(age~preBMIcat, data = baseline_anthro, weights = ~round(wt))

## For nulliparous
  # weighted count(percentage)
  sum(baseline_anthro$wt * baseline_anthro$Nulliparious) / sum(baseline_anthro$wt) * c(1, nrow(baseline_anthro))
  sum(baseline_anthro$wt[lBMI_idx] * baseline_anthro$Nulliparious[lBMI_idx]) / sum(baseline_anthro$wt[lBMI_idx]) * c(1, sum(baseline_anthro$wt[lBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  sum(baseline_anthro$wt[hBMI_idx] * baseline_anthro$Nulliparious[hBMI_idx]) / sum(baseline_anthro$wt[hBMI_idx]) * c(1, sum(baseline_anthro$wt[hBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  # test
  wtd.chi.sq(baseline_anthro$Nulliparious, baseline_anthro$preBMIcat, weight = baseline_anthro$wt)

## For Smoked
  # weighted count(percentage)
  sum(baseline_anthro$wt * baseline_anthro$Smoked) / sum(baseline_anthro$wt) * c(1, nrow(baseline_anthro))
  sum(baseline_anthro$wt[lBMI_idx] * baseline_anthro$Smoked[lBMI_idx]) / sum(baseline_anthro$wt[lBMI_idx]) * c(1, sum(baseline_anthro$wt[lBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  sum(baseline_anthro$wt[hBMI_idx] * baseline_anthro$Smoked[hBMI_idx]) / sum(baseline_anthro$wt[hBMI_idx]) * c(1, sum(baseline_anthro$wt[hBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  # invalid test with too small sample size
## For Alcohol
  #weighted count(percentage)
  sum(baseline_anthro$wt * (as.numeric(baseline_anthro$Alcohol_bin)-1)) / sum(baseline_anthro$wt) * c(1, nrow(baseline_anthro))
  sum(baseline_anthro$wt[lBMI_idx] * (as.numeric(baseline_anthro$Alcohol_bin)-1)[lBMI_idx]) / sum(baseline_anthro$wt[lBMI_idx]) * c(1, sum(baseline_anthro$wt[lBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  sum(baseline_anthro$wt[hBMI_idx] * (as.numeric(baseline_anthro$Alcohol_bin)-1)[hBMI_idx]) / sum(baseline_anthro$wt[hBMI_idx]) * c(1, sum(baseline_anthro$wt[hBMI_idx]) / sum(baseline_anthro$wt) * nrow(baseline_anthro))
  # test
  wtd.chi.sq(baseline_anthro$Alcohol_bin,baseline_anthro$preBMIcat,weight= baseline_anthro$wt)
  wtd.chi.sq(baseline_anthro$Alcohol,baseline_anthro$preBMIcat,weight= baseline_anthro$wt)
  

  ##### obtain the other anthro variables for all 321 selected individuals
lipid_anthro <- lipid %>%
  select(ID,age,GDM,wtPA12,race,preBMI,education,bw,bwz) %>%
  group_by(ID)

lipid_anthro <- lipid_anthro[!duplicated(lipid_anthro),]
lipid_anthro <- lipid_anthro %>%
  inner_join(select(anthro, -c("Age","preBMI")),by="ID")

lipid_anthro$preBMIcat <- as.factor(lipid_anthro$preBMI < 30)
lBMI_idx <- lipid_anthro$preBMI < 30
hBMI_idx <- lipid_anthro$preBMI >= 30
sum(lBMI_idx)
sum(hBMI_idx)
## For age
# median
  weighted.median(lipid_anthro$age,lipid_anthro$wtPA12)
  weighted.median(lipid_anthro$age[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx])
  weighted.median(lipid_anthro$age[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx])
# mean
  weighted.mean(lipid_anthro$age,lipid_anthro$wtPA12)
  sqrt(weighted.var(lipid_anthro$age,lipid_anthro$wtPA12))
# quantile
  weighted.quantile(lipid_anthro$age,lipid_anthro$wtPA12, probs = c(0.25,0.75))
  weighted.quantile(lipid_anthro$age[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx], probs = c(0.25,0.75))
  weighted.quantile(lipid_anthro$age[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx], probs = c(0.25,0.75))
# test
wilcox_test(age~preBMIcat, data = lipid_anthro, weights = ~round(wtPA12))

## For Race
  # test
  wtd.chi.sq(factor(lipid_anthro$race), factor(lipid_anthro$preBMIcat), weight = lipid_anthro$wtPA12)
  # weighted count and percentage
  for(k in 1:4){
    print(sum(lipid_anthro$wtPA12 * (as.numeric(factor(lipid_anthro$race == k)) - 1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro)))
  }
  for(k in 1:4){
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$race == k)) - 1)) / sum(lipid_anthro[lBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$race == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$race == k)) - 1)) / sum(lipid_anthro[hBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$race == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  }
## For Education
  # change education from 5 categories to 3 categories
  lipid_anthro$education <- ifelse(lipid_anthro$Education <= 2, 1, ifelse(lipid_anthro$Education >=4, 3,2))
  lipid_anthro$education <- factor(lipid_anthro$education, ordered = T)
  # test
  wtd.chi.sq(lipid_anthro$education, lipid_anthro$preBMIcat, weight = lipid_anthro$wtPA12)
  # weighted count and percentage
  for(k in 1:3){
    print(sum(lipid_anthro$wtPA12 * (as.numeric(factor(lipid_anthro$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro)))
  }
  for(k in 1:3){
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro[lBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro[hBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  }
  
  ## For nulliparous
  # weighted count(percentage)
  
  sum(lipid_anthro$wtPA12 * lipid_anthro$Nulliparious) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[lBMI_idx] * lipid_anthro$Nulliparious[lBMI_idx],na.rm = T) / 
    sum(lipid_anthro$wtPA12[lBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[lBMI_idx]) / 
                                             sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[hBMI_idx] * lipid_anthro$Nulliparious[hBMI_idx]) / 
    sum(lipid_anthro$wtPA12[hBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[hBMI_idx]) / 
                                             sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  # test
  wtd.chi.sq(lipid_anthro$Nulliparious, lipid_anthro$preBMIcat, weight = lipid_anthro$wtPA12)
  
  ## For Smoked
  # weighted count(percentage)
  sum(lipid_anthro$wtPA12 * lipid_anthro$Smoked) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[lBMI_idx] * lipid_anthro$Smoked[lBMI_idx]) / sum(lipid_anthro$wtPA12[lBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[lBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[hBMI_idx] * lipid_anthro$Smoked[hBMI_idx]) / sum(lipid_anthro$wtPA12[hBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[hBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  # invalid test with too small sample size
  ## For Alcohol
  #weighted count(percentage)
  sum(lipid_anthro$wtPA12 * (as.numeric(lipid_anthro$Alcohol_bin)-1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[lBMI_idx] * (as.numeric(lipid_anthro$Alcohol_bin)-1)[lBMI_idx]) / sum(lipid_anthro$wtPA12[lBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[lBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[hBMI_idx] * (as.numeric(lipid_anthro$Alcohol_bin)-1)[hBMI_idx]) / sum(lipid_anthro$wtPA12[hBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[hBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  # test
  wtd.chi.sq(lipid_anthro$Alcohol_bin,lipid_anthro$preBMIcat,weight= lipid_anthro$wtPA12)
  wtd.chi.sq(lipid_anthro$Alcohol,lipid_anthro$preBMIcat,weight= lipid_anthro$wtPA12)
  
  
## For preBMI
  # median
  weighted.median(lipid_anthro$preBMI,lipid_anthro$wtPA12)
  weighted.median(lipid_anthro$preBMI[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx])
  weighted.median(lipid_anthro$preBMI[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx])
  # mean and sd
  weighted.mean(lipid_anthro$preBMI,lipid_anthro$wtPA12)
  sqrt(weighted.var(lipid_anthro$preBMI,lipid_anthro$wtPA12))
  # quantile
  weighted.quantile(lipid_anthro$preBMI,lipid_anthro$wtPA12, probs = c(0.25,0.75))
  weighted.quantile(lipid_anthro$preBMI[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx], probs = c(0.25,0.75))
  weighted.quantile(lipid_anthro$preBMI[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx], probs = c(0.25,0.75))
  # test
  wilcox_test(preBMI~preBMIcat, data = lipid_anthro, weights = ~round(wtPA12))
  # recategorize preBMI into three groups 19 - 24.9, 25 - 29.9, and 30 - 45
  lipid_anthro$preBMIcat3 <- ifelse(lipid_anthro$preBMI <= 25, 1, ifelse(lipid_anthro$preBMI >=30,3,2))
  lipid_anthro$preBMIcat3 <- factor(lipid_anthro$preBMIcat3, ordered = T)
  # weighted count and percentage
  for(k in 1:3){
    print(sum(lipid_anthro$wtPA12 * (as.numeric(factor(lipid_anthro$preBMIcat3 == k)) - 1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro)))
  }
  for(k in 1:3){
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$preBMIcat3 == k)) - 1)) / sum(lipid_anthro[lBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$preBMIcat3 == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  }
  for(k in 1:3){
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$preBMIcat3 == k)) - 1)) / sum(lipid_anthro[hBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$preBMIcat3 == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  }
## For maternal GWG
  # median
  for(j in paste0("GWGv",c(0,1,2,4))){
    GWG = unlist(lipid_anthro[,j])
    print(weighted.median(GWG,lipid_anthro$wtPA12))
    print(weighted.median(GWG[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx]))
    print(weighted.median(GWG[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx]))
    # quantile
    print(weighted.quantile(GWG,lipid_anthro$wtPA12, probs = c(0.25,0.75)))
    print(weighted.quantile(GWG[lBMI_idx], lipid_anthro$wtPA12[lBMI_idx], probs = c(0.25,0.75)))
    print(weighted.quantile(GWG[hBMI_idx], lipid_anthro$wtPA12[hBMI_idx], probs = c(0.25,0.75)))
    # test
    print(wilcox_test(GWG~preBMIcat, data = lipid_anthro, weights = ~round(wtPA12)))
  }
  
# For GDM
  # weighted count(percentage)
  sum(lipid_anthro$wtPA12 * lipid_anthro$GDM) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[lBMI_idx] * lipid_anthro$GDM[lBMI_idx]) / sum(lipid_anthro$wtPA12[lBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[lBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  sum(lipid_anthro$wtPA12[hBMI_idx] * lipid_anthro$GDM[hBMI_idx]) / sum(lipid_anthro$wtPA12[hBMI_idx]) * c(1, sum(lipid_anthro$wtPA12[hBMI_idx]) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  wtd.chi.sq(lipid_anthro$GDM,lipid_anthro$preBMIcat,weight= lipid_anthro$wtPA12)
## For GDM, LGA, SGA. LBW, macrosomia
  # test
  wtd.chi.sq(lipid_anthro$education, lipid_anthro$preBMIcat, weight = lipid_anthro$wtPA12)
  # weighted count and percentage
  for(k in 1:3){
    print(sum(lipid_anthro$wtPA12 * (as.numeric(factor(lipid_anthro$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro)))
  }
  for(k in 1:3){
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro[lBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[lBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro[hBMI_idx,]$wtPA12))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * (as.numeric(factor(lipid_anthro[hBMI_idx,]$education == k)) - 1)) / sum(lipid_anthro$wtPA12) * nrow(lipid_anthro))
  }
  for(j in c("GDM","LGA","SGA","LBW","macro")){
    v <- unlist(lipid_anthro[,j])
    # test
    print(wtd.chi.sq(factor(v), lipid_anthro$preBMIcat, weight = lipid_anthro$wtPA12))
    # weighted count and percentage
    print(sum(lipid_anthro$wtPA12 * (as.numeric(factor(v))-1)) / sum(lipid_anthro$wtPA12) * c(1, nrow(lipid_anthro)))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * v[lBMI_idx],na.rm = T) / sum(lipid_anthro[lBMI_idx,]$wtPA12,na.rm = T))
    print(sum(lipid_anthro[lBMI_idx,]$wtPA12 * v[lBMI_idx],na.rm = T) / sum(lipid_anthro$wtPA12,na.rm = T) * nrow(lipid_anthro))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * v[hBMI_idx],na.rm = T) / sum(lipid_anthro[hBMI_idx,]$wtPA12,na.rm=T))
    print(sum(lipid_anthro[hBMI_idx,]$wtPA12 * v[hBMI_idx],na.rm = T) / sum(lipid_anthro$wtPA12,na.rm = T) * nrow(lipid_anthro))
  }

anthro_sel <- lipid_anthro %>%
  select(ID, age,GDM,race,preBMI,education,Education,bwz,bw.x,Length,Gender,SGA,LGA,PairID,SSF,sinceMealhr0,sinceMealhr1,sinceMealhr2,sinceMealhr4,ab_circ, head_circ,Nulliparious,LBW,macro,Smoked,Alcohol,Alcohol_bin,
         GWGv0,GWGv1,GWGv2,GWGv4,GAv0,GAv1,GAv2,GAv4,preBMIcat,preBMIcat3,measuredate)


## merge anthro info with the lipid info
lipid <- inner_join(lipid[,c(1,3,4,6:426,455)],anthro_sel)
lipid <- lipid[lipid$Smoked == 0,]
lipid$diet <- lipid$sinceMealhr0
lipid$diet[lipid$Visit==1] <- lipid$sinceMealhr1[lipid$Visit==1]
lipid$diet[lipid$Visit==2] <- lipid$sinceMealhr2[lipid$Visit==2]
lipid$diet[lipid$Visit==4] <- lipid$sinceMealhr4[lipid$Visit==4]
ggplot(lipid, aes(y = diet))+geom_histogram()+facet_wrap(~Visit)
aggregate(diet~Visit, data=lipid, function(x) quantile(x, c(0,0.05,0.25,0.5,0.75,0.95,1)))
write.csv(lipid,"Data/lipid_imputed_wosmoke.csv", row.names = F)
for(visit in c(0,1,2,4)){
  DF_InvN <- read.csv(paste0("tmp/lipid_InvN",visit,".csv"))
  DF_InvN <- inner_join(DF_InvN,anthro_sel,by="ID")
  DF_InvN <- DF_InvN[DF_InvN$Smoked==0,]
  write.csv(DF_InvN, paste0("tmp/lipid_InvN",visit,".csv"), row.names = F)
  DF_QN_InvN <- read.csv(paste0("tmp/lipid_QN_InvN",visit,".csv"))
  DF_QN_InvN <- inner_join(DF_QN_InvN,anthro_sel,by="ID")
  DF_QN_InvN <- DF_QN_InvN[DF_QN_InvN$Smoked==0,]
  write.csv(DF_QN_InvN,paste0("tmp/lipid_QN_InvN",visit,".csv"), row.names=F)
  
  DF_Std <- read.csv(paste0("tmp/lipid_Std",visit,".csv"))
  DF_Std <- inner_join(DF_Std,anthro_sel,by="ID")
  DF_Std <- DF_Std[DF_Std$Smoked==0,]
  write.csv(DF_Std, paste0("tmp/lipid_Std",visit,".csv"), row.names = F)
  DF_QN_Std <- read.csv(paste0("tmp/lipid_QN_Std",visit,".csv"))
  DF_QN_Std <- inner_join(DF_QN_Std,anthro_sel,by="ID")
  DF_QN_Std <- DF_Std[DF_QN_Std$Smoked==0,]
  write.csv(DF_QN_Std,paste0("tmp/lipid_QN_Std",visit,".csv"), row.names=F)
}

##----------------- Prepare data for trajectory analysis-------------------------

rm(list=ls())
comp.data4 <- read.csv("tmp/lipid_InvN4.csv")
for(i in 2:421){
  comp.data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all = T),
                      list(comp.data0[,c(1, 422, 424, 437, i)], 
                           comp.data1[,c(1, 422, 424, 437, i)], 
                           comp.data2[,c(1, 424, 437, i)], 
                           comp.data4[,c(1, 424, 437, i)]))
  comp.data <- data.frame(ID = comp.data$ID, PairID = coalesce(comp.data[,4], comp.data[,8], comp.data[,11], comp.data[,14]),
                          wt = coalesce(comp.data[,2], comp.data[,6]), V1 = comp.data[,5], V2 = comp.data[,9], V3 = comp.data[,12], V4 = comp.data[,15],
                          T1 = comp.data[,3], T2 = comp.data[,7], T3 = comp.data[,10], T4 = comp.data[,13])
  write.csv(comp.data, file = paste0("Data/sas_data/data",i-1,".csv"), row.names = F)
  
}
# for(i in 1:3){
#   comp.data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all = T),
#                       list(comp.data0[comp.data0$batch==i,c(1, 422, 424, 437, 4)], 
#                            comp.data1[comp.data1$batch==i,c(1, 422, 424, 437, 4)], 
#                            comp.data2[comp.data2$batch==i,c(1, 424, 437, 4)], 
#                            comp.data4[comp.data4$batch==i,c(1, 424, 437, 4)]))
#   comp.data <- data.frame(ID = comp.data$ID, PairID = coalesce(comp.data[,4], comp.data[,8], comp.data[,11], comp.data[,14]),
#                         wt = coalesce(comp.data[,2], comp.data[,6]), V1 = comp.data[,5], V2 = comp.data[,9], V3 = comp.data[,12], V4 = comp.data[,15],
#                         T1 = comp.data[,3], T2 = comp.data[,7], T3 = comp.data[,10], T4 = comp.data[,13])
#   write.csv(comp.data, file = paste0("Data/sas_data/data3_batch",i,".csv"), row.names = F)
# }

for(visit in c(0,1,2,4)){
  assign(paste0("lipid",visit), lipid[lipid$Visit == visit,])
}
for(i in 6:425){
  comp.data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all = T),
                      list(lipid0[,c(1, 455, 3, 431, 426, i)], 
                           lipid1[,c(1, 455, 3, 431, 426, i)], 
                           lipid2[,c(1, 3, 431, i)], 
                           lipid4[,c(1, 3, 431, i)]))
  comp.data <- data.frame(ID = comp.data$ID, PairID = coalesce(comp.data[,4], comp.data[,8], comp.data[,13], comp.data[,16]),
                        wt = coalesce(comp.data[,2], comp.data[,7]), batch = coalesce(comp.data[,5], comp.data[,10]),
                        V1 = comp.data[,6], V2 = comp.data[,11], V3 = comp.data[,14], V4 = comp.data[,17],
                        T1 = comp.data[,3], T2 = comp.data[,8], T3 = comp.data[,12], T4 = comp.data[,15])
  write.csv(comp.data, file = paste0("Data/sas_data/data", i - 5,".csv"), row.names = F)
}
# for(i in 1:3){
#   comp.data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all = T),
#                       list(lipid0[lipid0$batch==i,c(1, 455, 3, 431, 8)], 
#                            lipid1[lipid1$batch==i,c(1, 455, 3, 431, 8)], 
#                            lipid2[lipid2$batch==i,c(1, 3, 431, 8)], 
#                            lipid4[lipid4$batch==i,c(1, 3, 431, 8)]))
#   comp.data <- data.frame(ID = comp.data$ID, PairID = coalesce(comp.data[,4], comp.data[,8], comp.data[,11], comp.data[,14]),
#                           wt = coalesce(comp.data[,2], comp.data[,6]), V1 = comp.data[,5], V2 = comp.data[,9], V3 = comp.data[,12], V4 = comp.data[,15],
#                           T1 = comp.data[,3], T2 = comp.data[,7], T3 = comp.data[,10], T4 = comp.data[,13])
#   write.csv(comp.data, file = paste0("Data/sas_data/data3_ori_batch",i,".csv"), row.names = F)
# }
##---------------- Prepare anthro data for trajectory analysis------------------

rm(list = ls())
for(i in c(0,1,2,4)){
  DF <- read.csv(paste0("tmp/lipid_InvN",i,".csv"))
  assign(paste0("anthro",i),data.frame(ID = DF$ID, pairID = DF$PairID,wt = DF$wt,
                                       age = DF$age, race = factor(DF$race), education = factor(DF$Education, ordered = T), 
                                       preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = DF$Gender, SSF = DF$SSF,
                                       Circ = DF$Circ, Length = DF$Length, bwz = DF$bwz, bw = DF$bw.x/1000, alcohol = DF$Alcohol_bin, batch = factor(DF$batch),
                                       LGA = DF$LGA, SGA = DF$SGA, marco= DF$macro, LBW = DF$LBW))
}
anthro <- unique(rbind(anthro0[,-3],anthro1[,-3], anthro2[,-3], anthro4[,-3]))
anthro_wt <- unique(rbind(anthro0[,c(1,3)],anthro1[c(1,3)]))
anthro <- merge(anthro, anthro_wt, by = "ID", all.x = T)
write.csv(anthro, file = "traj_anthro.csv", row.names = F)
