rm(list = ls())
# setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
library(dplyr)
library(multcomp)
lipid <- read.csv("../Data/lipid_imputed_wosmoke_combat.csv", header = T)
anno <- read.csv("../Data/lipidomics_annotation_complete.csv",header=T, sep=",")
lipid_names <- data.frame(met_ID = substring(names(lipid[4:331]),4))
lipid_names <- inner_join(lipid_names,anno, by = "met_ID")
anthro <- read.csv("../traj_anthro.csv")
# anthro <- merge(anthro,unique(lipid[,c("ID","head_circ","ab_circ")]), by = "ID",
#                 all.x = T)
# anthro <- anthro[,-which(names(anthro) == "Circ")]
#write.csv(anthro,file = "traj_anthro.csv", row.names = F)
anthro$education <- factor(anthro$education > 3, ordered = T)
anthro$race <- factor(anthro$race == 3)
anthro$nulli <- factor(anthro$nulli)
anthro$alcohol <- factor(anthro$alcohol)
anthro$Gender <- factor(anthro$Gender)
##----------- Reorder the trajectory classes from most common to least----------##
factor_reorder <- function(x){
  return(factor(as.numeric(factor(x, levels = names(table(x))[order(table(x),decreasing = T)]))))
}

for(i in 1:328){
  out <- read.csv(paste0("Data/sas_data/res/out/out",i,".csv"))
  anthro <- merge(anthro,out[,c("ID","GROUP")], by = "ID")
  anthro$GROUP <- factor_reorder(anthro$GROUP)
  #anthro$GROUP <- factor(anthro$GROUP)
  names(anthro)[ncol(anthro)] <- paste0("id_",lipid_names$met_ID[i])
}

met_index <- 18:345
# summary the number of groups
num_groups <- rep(NA,328)
for(i in met_index){
  num_groups[i-17] <- length(levels(anthro[,i]))
}
table(num_groups)
# obtain the indices for 2-, 3-, 4-, 5-traj-class metabolites
for(k in 2:5){
  assign(paste0("traj_",k,"_idx"), which(num_groups==k))
}

#### --------- Summarize the results for 2-traj-class metabolites ------------------
traj_2_table <- data.frame(met_ID = names(anthro)[met_index[traj_2_idx]], 
                           beta0_1 = 0, beta1_1 = 0, beta2_1 = 0, betaCov_1 = 0,
                           beta0_2 = 0, beta1_2 = 0, beta2_2 = 0, betaCov_2 = 0)
ctt = 1
for(i in traj_2_idx){
  oe <- read.csv(paste0("Data/sas_data/res/oe/oe",i,".csv"))
  out <- read.csv(paste0("Data/sas_data/res/out/out",i,".csv"))
  most <- names(table(out$GROUP))[order(table(out$GROUP),decreasing = T)][1]
  if(most == "2"){
    traj_2_table$beta0_1[ctt] = oe[1, "INTERC02"]
    traj_2_table$beta0_2[ctt] = oe[1, "INTERC01"]
    if("LINEAR02" %in% names(oe)){
      traj_2_table$beta1_1[ctt] = oe[1, "LINEAR02"]
    }
    if("LINEAR01" %in% names(oe)){
      traj_2_table$beta1_2[ctt] = oe[1, "LINEAR01"]
    }
    if("QUADRA02" %in% names(oe)){
      traj_2_table$beta2_1[ctt] = oe[1, "QUADRA02"]
    }
    if("QUADRA01" %in% names(oe)){
      traj_2_table$beta2_2[ctt] = oe[1, "QUADRA01"]
    }
    
    traj_2_table$betaCov_1[ctt] = oe[1,"TCOV1_NUm02"]
    traj_2_table$betaCov_2[ctt] = oe[1,"TCOV1_NUm01"]
  }else{
    traj_2_table$beta0_1[ctt] = oe[1, "INTERC01"]
    traj_2_table$beta0_2[ctt] = oe[1, "INTERC02"]
    if("LINEAR01" %in% names(oe)){
      traj_2_table$beta1_1[ctt] = oe[1, "LINEAR01"]
    }
    if("LINEAR02" %in% names(oe)){
      traj_2_table$beta1_2[ctt] = oe[1, "LINEAR02"]
    }
    if("QUADRA01" %in% names(oe)){
      traj_2_table$beta2_1[ctt] = oe[1, "QUADRA01"]
    }
    if("QUADRA02" %in% names(oe)){
      traj_2_table$beta2_2[ctt] = oe[1, "QUADRA02"]
    }
    
    traj_2_table$betaCov_1[ctt] = oe[1,"TCOV1_NUm01"]
    traj_2_table$betaCov_2[ctt] = oe[1,"TCOV1_NUm02"]
  }
  
  ctt = ctt + 1
}
## add the 2 traj-class columns indicating the trend of the trajectory
get_trajClass <- function(x){
  # Find the trend of the curve in the range from 0-40 (whole period of pregnancy)
    beta1 <- x[1]
    beta2 <- x[2]
    det <- - beta1 / (2*beta2)
    if(beta2 == 0){
      if(beta1 == 0){
        return("Constant")
      }else if(beta1 < 0){
        return("Linear decreasing")
      }else{
        return("Linear increasing")
      }
    }else if(beta2 > 0){
      if(det < 0){
        return("Quadratic increasing")
      }else if(det > 40){
        return("Quadratic decreasing")
      }else{
        return("Decreasing - increasing")
      }
    }else{
      if(det < 0){
        return("Quadratic decreasing")
      }else if(det > 40){
        return("Quadratic increasing")
      }else{
        return("Increasing - decreasing")
      }
    }
    
}
add_traj <- function(table, Class){
  ctt <- 3
  for(k in 1:Class){
    table <- cbind(table, factor(apply(table[,c(ctt,ctt+1)],1,get_trajClass),
                                 levels = c("Constant","Linear increasing","Quadratic increasing",
                                            "Linear decreasing","Quadratic decreasing",
                                            "Increasing - decreasing","Decreasing - increasing")))
    names(table)[ncol(table)] <- paste0("traj",k)
    ctt = ctt + 4
  }
  return(table)
}
traj_2_table$traj1 <- factor(apply(traj_2_table[,c("beta1_1","beta2_1")],1,get_trajClass),
                             levels = c("Constant","Linear increasing","Quadratic increasing",
                                        "Linear decreasing","Quadratic decreasing",
                                        "Increasing - decreasing","Decreasing - increasing"))
traj_2_table$traj2 <- factor(apply(traj_2_table[,c("beta1_2","beta2_2")],1,get_trajClass),
                             levels = c("Constant","Linear increasing","Quadratic increasing",
                                        "Linear decreasing","Quadratic decreasing",
                                        "Increasing - decreasing","Decreasing - increasing"))
library(lmerTest)
library(car)
##----------- GLM of anthro outcome vs. trajectory class-------------##
#Define LME function: random intercept - pairID, weighted
anthros <- c("bw", "bwz", "SSF", "Length","head_circ","ab_circ")
anthros_idx <- which(names(anthro) %in% anthros)
anthros <- names(anthro)[anthros_idx]
ctt = 1
for(idx in anthros_idx){
  DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                   age = anthro$age, race = anthro$race, education = anthro$education, 
                   preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                   outcome = anthro[,idx], anthro[,met_index])
  lme_met <- function(x) {
    nDF2 <- DF[,c(1:10)]
    nDF2$met <- x
    nDF2 <- na.omit(nDF2)
    #model = lme(anthro ~ met + age + as.factor(race) + preBMI + parity + female + GADel + glucose, 
    #            data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    model = lmer(outcome ~ met + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
                 data = nDF2, weights = wt)
    est <- c(summary(model)$coefficients[2,], Anova(model,type = 2)[[3]][1])
    return(est)
  }
  res <- as.data.frame(t(sapply(anthro[,met_index[traj_2_idx]], lme_met)))[,c(1,2,5,6)]
  colnames(res) <- c("beta","se","pval_t","pval_2")
  res$met_ID <- row.names(res)
  res <- merge(res, traj_2_table, by = "met_ID", all.x=T)
  write.csv(res,paste0("tmp/traj_",anthros[ctt],"_2class.csv"),row.names = F)
  ctt = ctt + 1
}

#### --------- Summarize the results for 3-traj-class metabolites ------------------
traj_3_table <- data.frame(met_ID = names(anthro)[met_index[traj_3_idx]], 
                           beta0_1 = 0, beta1_1 = 0, beta2_1 = 0, betaCov_1 = 0,
                           beta0_2 = 0, beta1_2 = 0, beta2_2 = 0, betaCov_2 = 0,
                           beta0_3 = 0, beta1_3 = 0, beta2_3 = 0, betaCov_3 = 0)

get_estimation <- function(oe, k){
  beta <- rep(0,4)
  beta[1] <- oe[1,paste0("INTERC0",k)]
  if(paste0("LINEAR0",k)%in% names(oe)){beta[2] <- oe[1,paste0("LINEAR0",k)]}
  if(paste0("QUADRA0",k)%in% names(oe)){beta[3] <- oe[1,paste0("QUADRA0",k)]}
  beta[4] <- oe[1,paste0("TCOV1_NUm0",k)]
  return(beta)
}


ctt = 1
for(i in traj_3_idx){
  oe <- read.csv(paste0("Data/sas_data/res/oe/oe",i,".csv"))
  out <- read.csv(paste0("Data/sas_data/res/out/out",i,".csv"))
  factor_order <- names(table(out$GROUP))[order(table(out$GROUP),decreasing = T)]
  col_stat <- 2
  for(k in factor_order){
    traj_3_table[ctt,col_stat:(col_stat + 3)] <- get_estimation(oe,k)
    col_stat <- col_stat + 4
  }
  ctt = ctt + 1
}
## add the 3 traj-class columns indicating the trend of the trajectory
traj_3_table <- add_traj(traj_3_table,3)

##----------- GLM of anthro outcome vs. trajectory class-------------##
#Define LME function: random intercept - pairID, weighted
anthros <- c("bw", "bwz", "SSF", "Length","head_circ","ab_circ")
anthros_idx <- which(names(anthro) %in% anthros)
anthros <- names(anthro)[anthros_idx]
ctt = 1
for(idx in anthros_idx){
  DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                   age = anthro$age, race = anthro$race, education = anthro$education, 
                   preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                   outcome = anthro[,idx], anthro[,met_index])
  lme_met <- function(x) {
    nDF2 <- DF[,c(1:10)]
    nDF2$met <- x
    nDF2 <- na.omit(nDF2)
    #model = lme(anthro ~ met + age + as.factor(race) + preBMI + parity + female + GADel + glucose, 
    #            data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    model = lmer(outcome ~ met + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
                 data = nDF2, weights = wt)
    postHocs <- glht(model, linfct = mcp("met" = "Tukey"))
    est <- c(Anova(model,type = 2)[[3]][1], summary(postHocs)$test$coefficients,
             summary(postHocs)$test$sigma,
             summary(postHocs)$test$pvalues[1:3])
    
    return(est)
  }
  res <- as.data.frame(t(sapply(anthro[,met_index[traj_3_idx]], lme_met)))
  names(res)[c(1,5:10)] <- c("pval_2",paste0("se_comp",1:3), paste0("pval_comp",1:3))
  res$met_ID <- row.names(res)
  res <- merge(res, traj_3_table, by = "met_ID", all.x=T)
  write.csv(res,paste0("tmp/traj_",anthros[ctt],"_3class.csv"),row.names = F)
  ctt = ctt + 1
}

#### --------- Summarize the results for 4-traj-class metabolites ------------------
traj_4_table <- data.frame(met_ID = names(anthro)[met_index[traj_4_idx]], 
                           beta0_1 = 0, beta1_1 = 0, beta2_1 = 0, betaCov_1 = 0,
                           beta0_2 = 0, beta1_2 = 0, beta2_2 = 0, betaCov_2 = 0,
                           beta0_3 = 0, beta1_3 = 0, beta2_3 = 0, betaCov_3 = 0,
                           beta0_4 = 0, beta1_4 = 0, beta2_4 = 0, betaCov_4 = 0)
ctt = 1
for(i in traj_4_idx){
  oe <- read.csv(paste0("Data/sas_data/res/oe/oe",i,".csv"))
  out <- read.csv(paste0("Data/sas_data/res/out/out",i,".csv"))
  factor_order <- names(table(out$GROUP))[order(table(out$GROUP),decreasing = T)]
  col_stat <- 2
  for(k in factor_order){
    traj_4_table[ctt,col_stat:(col_stat + 3)] <- get_estimation(oe,k)
    col_stat <- col_stat + 4
  }
  ctt = ctt + 1
}
## add the 4 traj-class columns indicating the trend of the trajectory
traj_4_table <- add_traj(traj_4_table, 4)
##----------- GLM of anthro outcome vs. trajectory class-------------##
#Define LME function: random intercept - pairID, weighted
anthros_idx <- which(names(anthro) %in% anthros)
anthros <- names(anthro)[anthros_idx]
ctt = 1
for(idx in anthros_idx){
  DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                   age = anthro$age, race = anthro$race, education = anthro$education, 
                   preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                   outcome = anthro[,idx], anthro[,met_index])
  lme_met <- function(x) {
    nDF2 <- DF[,c(1:10)]
    nDF2$met <- x
    nDF2 <- na.omit(nDF2)
    #model = lme(anthro ~ met + age + as.factor(race) + preBMI + parity + female + GADel + glucose, 
    #            data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    model = lmer(outcome ~ met + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
                 data = nDF2, weights = wt)
    postHocs <- glht(model, linfct = mcp("met" = "Tukey"))
    est <- c(Anova(model,type = 2)[[3]][1], summary(postHocs)$test$coefficients,
             summary(postHocs)$test$sigma,
             summary(postHocs)$test$pvalues[1:6])
    
    return(est)
  }
  res <- as.data.frame(t(sapply(anthro[,met_index[traj_4_idx]], lme_met)))
  names(res)[c(1,8:19)] <- c("pval_2",paste0("se_comp",1:6), paste0("pval_comp",1:6))
  res$met_ID <- row.names(res)
  res <- merge(res, traj_4_table, by = "met_ID", all.x=T)
  write.csv(res,paste0("tmp/traj_",anthros[ctt],"_4class.csv"),row.names = F)
  ctt = ctt + 1
}

#### --------- Summarize the results for 5-traj-class metabolites ------------------
traj_5_table <- data.frame(met_ID = names(anthro)[met_index[traj_5_idx]], 
                           beta0_1 = 0, beta1_1 = 0, beta2_1 = 0, betaCov_1 = 0,
                           beta0_2 = 0, beta1_2 = 0, beta2_2 = 0, betaCov_2 = 0,
                           beta0_3 = 0, beta1_3 = 0, beta2_3 = 0, betaCov_3 = 0,
                           beta0_4 = 0, beta1_4 = 0, beta2_4 = 0, betaCov_4 = 0,
                           beta0_5 = 0, beta1_5 = 0, beta2_5 = 0, betaCov_5 = 0)
ctt = 1
for(i in traj_5_idx){
  oe <- read.csv(paste0("Data/sas_data/res/oe/oe",i,".csv"))
  out <- read.csv(paste0("Data/sas_data/res/out/out",i,".csv"))
  factor_order <- names(table(out$GROUP))[order(table(out$GROUP),decreasing = T)]
  col_stat <- 2
  for(k in factor_order){
    traj_5_table[ctt,col_stat:(col_stat + 3)] <- get_estimation(oe,k)
    col_stat <- col_stat + 4
  }
  ctt = ctt + 1
}
## add the 5 traj-class columns indicating the trend of the trajectory
traj_5_table <- add_traj(traj_5_table, 5)
##----------- GLM of anthro outcome vs. trajectory class-------------##
#Define LME function: random intercept - pairID, weighted
for(ctt in 1:6){
  idx = anthros_idx[ctt]
  DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                   age = anthro$age, race = anthro$race, education = anthro$education, 
                   preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                   outcome = anthro[,idx], anthro[,met_index])
  lme_met <- function(x) {
    nDF2 <- DF[,c(1:10)]
    nDF2$met <- x
    nDF2 <- na.omit(nDF2)
    #model = lme(anthro ~ met + age + as.factor(race) + preBMI + parity + female + GADel + glucose, 
    #            data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    model = lmer(outcome ~ met + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
                 data = nDF2, weights = wt)
    postHocs <- glht(model, linfct = mcp("met" = "Tukey"))
    est <- c(Anova(model,type = 2)[[3]][1], summary(postHocs)$test$coefficients,
             summary(postHocs)$test$sigma,
             summary(postHocs)$test$pvalues[1:10])
    
    return(est)
  }
  res <- as.data.frame(t(sapply(anthro[,met_index[traj_5_idx]], lme_met)))
  names(res)[c(1,12:31)] <- c("pval_2",paste0("se_comp",1:10), paste0("pval_comp",1:10))
  res$met_ID <- row.names(res)
  res <- merge(res, traj_5_table, by = "met_ID", all.x=T)
  write.csv(res,paste0("tmp/traj_",anthros[ctt],"_5class.csv"),row.names = F)
}


#### --------- Merge F test results of all trajectory patterns together ------
F_pval <- NULL
outcome_list = c("Sum of skinfolds", "Length", "Birthweight Z-score",
                 "Birth weight","Head circumference","Abdominal circumference")
for(idx in 1:6){
  ant = anthros[idx]
  res <- vector(mode="list", length = 4)
  pval <- NULL
  for(k in 2:5){
    res_tmp <- read.csv(paste0("tmp/traj_",ant,"_",k,"class.csv"))
    pval <- c(pval,res_tmp$pval_2)
    res[[k-1]] <- res_tmp
  }
  pFDR <- p.adjust(pval, method = "fdr")
  pBonf <- p.adjust(pval, method = "bonferroni")
  n_class_start <- cumsum(c(1,lapply(res,nrow)))
  n_class_end <- cumsum(c(lapply(res,nrow)))
  for(i in 1:4){
    res[[i]]$pFDR <- pFDR[n_class_start[i] : n_class_end[i]]
    res[[i]]$pBonf <- pBonf[n_class_start[i] : n_class_end[i]]
    F_pval = rbind(F_pval, data.frame(met_ID = substring(res[[i]]$met_ID,4), pval = res[[i]]$pval_2,
                                      pFDR = res[[i]]$pFDR, pBonf = res[[i]]$pBonf, outcome = outcome_list[idx]))
    write.csv(res[[i]],file = paste0("tmp/traj_",ant,"_",i + 1,"class_adjp.csv"),
              row.names = F)
  }
}
nrow(F_pval)
F_pval = inner_join(F_pval,anno,by="met_ID")
F_pval$significance <- F_pval$pFDR < 0.05
F_pval$label <- ifelse(F_pval$significance, F_pval$Annotation, NA)

   
  anno.tab <- unique(F_pval[,c("Annotation","Lipid.class")])
  anno.tab <- table(anno.tab$Annotation, anno.tab$Lipid.class)
  anno.min <- (1:nrow(anno.tab))-0.5
  anno.max <- (1:nrow(anno.tab))+0.5
  min.list <- apply(anno.tab, 2, function(x) min(anno.min[x > 0]))
  max.list <- apply(anno.tab, 2, function(x) max(anno.max[x > 0]))
  x.anno <- data.frame(fill = as.factor(names(min.list)), xmin = min.list, xmax = max.list)
  ggplot(F_pval, aes(x = Annotation, y = -log10(pFDR))) + 
    geom_point(aes(col = outcome,shape = outcome,alpha = as.numeric(significance)-1.6),size = 3)+
    geom_hline(yintercept = -log10(0.05)) +
    geom_rect(x.anno, mapping = aes(xmin = xmin, xmax = xmax, fill=as.factor(fill)), 
              ymin = -Inf, ymax=Inf, alpha = 0.15, inherit.aes = F) +
    scale_fill_manual(name="Lipid Class",
                      limits = x.anno$fill,
                      labels = x.anno$fill,
                      values= mycolors12)+
    scale_color_discrete(name = "Fetal outcome")+
    scale_shape_discrete(name = "Fetal outcome")+
    scale_alpha_continuous(guide = F)+
    scale_x_discrete(name="Metabolite annotation",breaks=NULL)+
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 10)+theme_bw()

  # plot outcomes as strata and sperate them into different panels
  F_pval$label_comp = ifelse(F_pval$significance,F_pval$Annotation_comp,NA)
  anno.tab <- unique(F_pval[,c("Annotation","Lipid.class")])
  anno.tab <- table(anno.tab$Annotation, anno.tab$Lipid.class)
  anno.min <- (1:nrow(anno.tab))-0.5
  anno.max <- (1:nrow(anno.tab))+0.5
  min.list <- apply(anno.tab, 2, function(x) min(anno.min[x > 0]))
  max.list <- apply(anno.tab, 2, function(x) max(anno.max[x > 0]))
  x.anno <- data.frame(fill = as.factor(names(min.list)), xmin = min.list, xmax = max.list)
  ggplot(F_pval, aes(x = Annotation, y = -log10(pFDR))) + 
    geom_point(aes(col = Lipid.class,fill = Lipid.class, alpha = as.numeric(significance)-1.6),size = 3)+
    geom_hline(yintercept = -log10(0.05)) +
    geom_rect(x.anno, mapping = aes(xmin = xmin, xmax = xmax, fill=as.factor(fill)), 
              ymin = -Inf, ymax=Inf, alpha = 0.15, inherit.aes = F) +
    scale_alpha_continuous(guide = F)+
    scale_x_discrete(name="Metabolite annotation",breaks=NULL)+
    geom_text_repel(aes(label = label),size = 3,
                    max.overlaps = 5,
                    min.segment.length = Inf,
                    )+
    theme_bw()+
    facet_wrap(~ factor(outcome))+
    scale_fill_manual(name = "Lipid Class",
                       values= mycolors12) +
    scale_color_manual(name="Lipid Class",
                      values= mycolors12)
    
    F_pval_sig <- F_pval[F_pval$significance == T,]
    write.csv(F_pval_sig,file="Figures/draft/trajectory analysis/Ftest_pval.csv", row.names = F)
    
    
#### --------- Plot for different patterns ---------------
    ## ------------- Two-class -----------------------
    k = 2
    plot.data=NULL
    for(i in 1:6){
      ant = anthros[i]
      data <- read.csv(paste0("tmp/traj_",ant,"_",k,"class_adjp.csv"))
      data$outcome = outcome_list[i]
      plot.data <- rbind(plot.data,data)
    }
    plot.data$met_ID = substring(plot.data$met_ID,4)
    plot.data <- inner_join(plot.data, anno, by="met_ID")
    plot.data$significance <- plot.data$pFDR < 0.05
    plot.data$label <- ifelse(plot.data$significance,plot.data$Annotation,NA)
    plot.data$traj1_overall <- plot.data$traj1
    plot.data$traj1_overall[plot.data$traj1_overall %in% c("Linear increasing", "Quadratic increasing")] <- "Increasing"
    plot.data$traj1_overall[plot.data$traj1_overall %in% c("Linear decreasing", "Quadratic decreasing")] <- "Decreasing"
    plot.data$traj2_overall <- plot.data$traj2
    plot.data$traj2_overall[plot.data$traj2_overall %in% c("Linear increasing", "Quadratic increasing")] <- "Increasing"
    plot.data$traj2_overall[plot.data$traj2_overall %in% c("Linear decreasing", "Quadratic decreasing")] <- "Decreasing"
    
    ggplot(plot.data, aes(x = beta, y = -log10(pFDR))) + 
      geom_point(aes(col = Lipid.class, shape= outcome,alpha = as.numeric(significance) - 1.4), size = 2.5)+
      geom_text_repel(aes(label = label),size = 2)+
      geom_vline(xintercept = 0, col = "grey")+
      geom_hline(yintercept = -log10(0.05),col="grey")+
      scale_alpha_continuous(guide = F) + 
      scale_color_manual(values = mycolors12,
                           name = "Lipid class")+
      scale_shape_discrete(name = "Fetal outcome")+
      facet_grid(traj2_overall~traj1_overall)+
      theme_classic()
    write.csv(plot.data[,c("outcome","met_ID","beta","se",
                       "pval_2","traj1","traj2","traj1_overall","traj2_overall",
                       "pFDR","pBonf","significance","Annotation","Lipid.class",
                       "label")], file = "Figures/draft/trajectory analysis/multcomp_2class.csv", row.names = F)
    
    ## ------------- Three-class ---------------
    k = 3
    plot.data=NULL
    for(i in 1:6){
      ant = anthros[i]
      data <- read.csv(paste0("tmp/traj_",ant,"_",k,"class_adjp.csv"))
      data$outcome = outcome_list[i]
      plot.data <- rbind(plot.data,data)
    }
    names(plot.data)[3:5] <- paste0("beta_comp",1:3)
    plot.data$met_ID = substring(plot.data$met_ID,4)
    plot.data <- inner_join(plot.data, anno, by="met_ID")
   
    plot.data$comparison1 <- paste0(plot.data$traj2, " \n vs. \n ", plot.data$traj1)
    plot.data$comparison2 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj1)
    plot.data$comparison3 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj2)
    plot.data$overall_sig <- plot.data$pFDR < 0.05
    
    plot <- data.frame (beta = unlist(plot.data[,3:5]),
                        se = unlist(plot.data[,6:8]),
                        pval = unlist(plot.data[,9:11]),
                        comparison = as.factor(unlist(plot.data[,42:44])),
                        pFDR = rep(plot.data$pFDR,3),
                        Annotation = rep(plot.data$Annotation,3),
                        Lipid.class = rep(plot.data$Lipid.class,3),
                        outcome = rep(plot.data$outcome,3)
                        )
    plot$significance <- (plot$pFDR < 0.05) & (plot$pval < 0.05)
    plot$label <- ifelse(plot$significance,plot$Annotation,NA)
    write.csv(plot, file = "Figures/draft/trajectory analysis/multcomp_3class.csv",
              row.names = F)
    #sigPerComp <- apply(table(plot$significance, plot$comparison, plot$outcome),c(1,2),
    #       function(x) max(x))[2,]
    #sigPerComp <- data.frame(comparison = names(sigPerComp), number = as.vector(sigPerComp))
    #anno_break_sig <- unique(plot$Annotation[plot$significance])
   p <- ggplot(subset(plot,significance), aes(x = beta, y = Annotation))+
      geom_pointrange(aes(col = Lipid.class,xmin = beta-se, xmax = beta+se, alpha = as.numeric(significance) - 1.4))+
      geom_text_repel(aes(label = label), size = 3) + 
      scale_alpha_continuous(guide = F)+
      theme_classic() +
      geom_vline(xintercept = 0) +
      scale_color_manual(values= mycolors12,
                         name = "Lipid class")+
      scale_y_discrete(aes(limits = label),
                       guide = guide_axis(check.overlap = T))+
      facet_grid(as.factor(comparison)~outcome,
                 scales = "free")+
      theme(legend.position = "bottom",
            strip.text.y.right = element_text(size = 8,angle = 0),
            axis.title.y = element_blank())
  print(p)  
  
    ## ------------- Four-class ------------
  k = 4
  plot.data=NULL
  for(i in 1:6){
    ant = anthros[i]
    data <- read.csv(paste0("tmp/traj_",ant,"_",k,"class_adjp.csv"))
    data$outcome = outcome_list[i]
    plot.data <- rbind(plot.data,data)
  }
  names(plot.data)[3:8] <- paste0("beta_comp",1:6)
  plot.data$met_ID = substring(plot.data$met_ID,4)
  plot.data <- inner_join(plot.data, anno, by="met_ID")
  
  plot.data$comparison1 <- paste0(plot.data$traj2, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison2 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison3 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison4 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj2)
  plot.data$comparison5 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj2)
  plot.data$comparison6 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj3)
  plot.data$overall_sig <- plot.data$pFDR < 0.05
  
  plot <- data.frame (beta = unlist(plot.data[,3:8]),
                      se = unlist(plot.data[,9:14]),
                      pval = unlist(plot.data[,15:20]),
                      comparison = as.factor(unlist(plot.data[,56:61])),
                      pFDR = rep(plot.data$pFDR,6),
                      Annotation = rep(plot.data$Annotation,6),
                      Lipid.class = rep(plot.data$Lipid.class,6),
                      outcome = rep(plot.data$outcome,6)
  )
  plot$significance <- (plot$pFDR < 0.05) & (plot$pval < 0.05)
  plot$label <- ifelse(plot$significance,plot$Annotation,NA)
  write.csv(plot, file = "Figures/draft/trajectory analysis/multcomp_4class.csv",
            row.names = F)
  #sigPerComp <- apply(table(plot$significance, plot$comparison, plot$outcome),c(1,2),
  #       function(x) max(x))[2,]
  #sigPerComp <- data.frame(comparison = names(sigPerComp), number = as.vector(sigPerComp))
  #anno_break_sig <- unique(plot$Annotation[plot$significance])
  p <- ggplot(subset(plot,significance), aes(x = beta, y = Annotation))+
    geom_pointrange(aes(col = Lipid.class,xmin = beta-se, xmax = beta+se, alpha = as.numeric(significance) - 1.4))+
    geom_text_repel(aes(label = label), size = 3) + 
    scale_alpha_continuous(guide = F)+
    theme_classic() +
    geom_vline(xintercept = 0) +
    scale_color_manual(values= mycolors12,
                       name = "Lipid class")+
    scale_y_discrete(aes(limits = label),
                     guide = guide_axis(check.overlap = T))+
    facet_grid(as.factor(comparison)~outcome,
               scales = "free")+
    theme(legend.position = "bottom",
          strip.text.y.right = element_text(size = 8,angle = 0),
          axis.title.y = element_blank())
  print(p)  
  
    ## ------------- Five-class ------------
  k = 5
  plot.data=NULL
  for(i in 1:6){
    ant = anthros[i]
    data <- read.csv(paste0("tmp/traj_",ant,"_",k,"class_adjp.csv"))
    data$outcome = outcome_list[i]
    plot.data <- rbind(plot.data,data)
  }
  names(plot.data)[3:12] <- paste0("beta_comp",1:10)
  plot.data$met_ID = substring(plot.data$met_ID,4)
  plot.data <- inner_join(plot.data, anno, by="met_ID")
  
  plot.data$comparison1 <- paste0(plot.data$traj2, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison2 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison3 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison4 <- paste0(plot.data$traj5, " \n vs. \n ", plot.data$traj1)
  plot.data$comparison5 <- paste0(plot.data$traj3, " \n vs. \n ", plot.data$traj2)
  plot.data$comparison6 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj2)
  plot.data$comparison7 <- paste0(plot.data$traj5, " \n vs. \n ", plot.data$traj2)
  plot.data$comparison8 <- paste0(plot.data$traj4, " \n vs. \n ", plot.data$traj3)
  plot.data$comparison9 <- paste0(plot.data$traj5, " \n vs. \n ", plot.data$traj3)
  plot.data$comparison10 <- paste0(plot.data$traj5, " \n vs. \n ", plot.data$traj4)
  plot.data$overall_sig <- plot.data$pFDR < 0.05
  
  plot <- data.frame (beta = unlist(plot.data[,3:12]),
                      se = unlist(plot.data[,13:22]),
                      pval = unlist(plot.data[,23:32]),
                      comparison = as.factor(unlist(plot.data[,73:82])),
                      pFDR = rep(plot.data$pFDR,10),
                      Annotation = rep(plot.data$Annotation,10),
                      Lipid.class = rep(plot.data$Lipid.class,10),
                      outcome = rep(plot.data$outcome,10)
  )
  plot$significance <- (plot$pFDR < 0.05) & (plot$pval < 0.05)
  plot$label <- ifelse(plot$significance,plot$Annotation,NA)
  write.csv(plot, file = "Figures/draft/trajectory analysis/multcomp_5class.csv",
            row.names = F)
  #sigPerComp <- apply(table(plot$significance, plot$comparison, plot$outcome),c(1,2),
  #       function(x) max(x))[2,]
  #sigPerComp <- data.frame(comparison = names(sigPerComp), number = as.vector(sigPerComp))
  #anno_break_sig <- unique(plot$Annotation[plot$significance])
  p <- ggplot(subset(plot,significance), aes(x = beta, y = Annotation))+
    geom_pointrange(aes(col = Lipid.class,xmin = beta-se, xmax = beta+se, alpha = as.numeric(significance) - 1.4))+
    geom_text_repel(aes(label = label), size = 3) + 
    scale_alpha_continuous(guide = F)+
    theme_classic() +
    geom_vline(xintercept = 0) +
    scale_color_manual(values= mycolors12,
                       name = "Lipid class")+
    scale_y_discrete(aes(limits = label),
                     guide = guide_axis(check.overlap = T))+
    facet_grid(as.factor(comparison)~outcome,
               scales = "free")+
    theme(legend.position = "bottom",
          strip.text.y.right = element_text(size = 8,angle = 0),
          axis.title.y = element_blank())
  print(p)  
  