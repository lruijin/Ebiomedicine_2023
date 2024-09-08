setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
rm(list= ls())
lipid <- read.csv("../Data/lipid_imputed_wosmoke_combat.csv", header = T)
anno <- read.csv("../Data/lipidomics_annotation_complete.csv",header=T, sep=",")
library(lme4)
library(lmerTest)
library(dplyr)
met_idx = 4:423
lipid$preBMIcat <- as.numeric(lipid$preBMIcat)

anthros <- c("bw","bwz","Length","ab_circ","head_circ","SSF")
for (visit in c(0, 1, 2, 4)) {    # By visit
  for(j in 1:length(anthros)){
    ###1.Read and prepare data---------------------------------------------------------------
    DF <- lipid[lipid$Visit == visit & lipid$preBMI >= 30,]
    if(j == 1){
      DF <- data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                       age = DF$age, race = factor(DF$race==3), education = factor(DF$Education > 3, ordered = T), 
                       preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                       anthro = DF$bw.x/1000, 
                       DF[,met_idx])
    }else{
      anthro_idx = which(names(DF) == anthros[j])
      DF <- data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                       age = DF$age, race = factor(DF$race==3), education = factor(DF$Education > 3, ordered = T), 
                       preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                       anthro = DF[,anthro_idx],
                       DF[,met_idx])
    }
    #Cor <- cor(DF[,15:342],DF[,c("anthro","age","race","education","GA_blood","Gender","nulli","alcohol","post")],
    #           weights.x = 1/DF$v, weights.y = 1/DF$v, use="p")
    #pTable <- corPvalueFisher(Cor, nrow(DF))
    #apply(pTable,2,summary)
    # met_index_df <- c(15:434) # all lipids
    met_index_df <- c(12:339) # all known lipids 328 intotal
    
    ###2.Estimate batch-specific effects-----------------------------------------------------
    #Define LME function: random intercept - pairID, weighted
    lme_met <- function(x) {
      DF2 <- DF[,c(1:11)]
      DF2$met <- x
      DF2 <- na.omit(DF2)
      model = lmer(anthro ~ met + age + race + education + GA_blood + Gender + nulli + (1|pairID) , 
                  data = DF2, weights = 1/v, control = lmerControl(optimizer = "optimx", optCtrl = list(method="nlminb")))
      #model = lme(anthro ~ met + age + race + education + GA_blood + Gender + nulli, 
      #            data = DF2, random = ~1 | pairID, weights = varFixed(~v))
      if(any( grepl("failed to converge", model@optinfo$conv$lme4$messages) )){
        est = c(summary(model)$coefficients[2,],F)
      }else{
        est <- c(summary(model)$coefficients[2,],T)
      }
      return(est)
    }
    
    
    res <- as.data.frame(t(sapply(DF[,met_index_df], lme_met)))[,c(1,2,5,6)]
    colnames(res) <- c("beta","se","pval","convergence")
    res$met_ID <- substring(row.names(res), 4)
    res <- inner_join(anno, res, by = "met_ID")
    res$pBonf <- p.adjust(res$pval,"bonferroni")
    res$pFDR <- p.adjust(res$pval,"fdr")
    
    write.csv(res, paste("tmp/known/", anthros[j], visit, "_sBMI.csv",sep=""), row.names=F)
  }
}

## summary the result:
for(i in c(0,1,2,4)){
  assign(paste0("res",i), read.csv(paste0("tmp/known/head_circ",i,"_lBMI",".csv"), header=T, sep=","))
  assign(paste0("res",i), data.frame(get(paste0("res",i)), Visit = i))
}

res = rbind(res0,res1,res2,res4)
sum(res$convergence==0)

write.csv(res,"Result/Table 2/individual analysis_continuous/head_circ/table2_head_circ_long_lBMI.csv",row.names = F)

tab.res = res %>% filter(pFDR < 0.05) %>%
  arrange(Annotation, met_ID) %>%
  select(Anno.,Annotation,met_ID,Visit,beta,se,pval,pFDR)

res_merge <- merge(merge(merge(res0[,c("Annotation", "met_ID","beta","se","pval","pBonf","pFDR")],
                               res1[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID", all.x = T),
                               res2[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID", all.x = T),
                               res4[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID", all.x = T)

names(res_merge)[c(3:22)] <- paste0(rep(c("beta","se","pval","pBonf","pFDR"),4),rep(c(0,1,2,4),each=5))
res_merge$Ann. <- substring(res_merge$Annotation,1,3)
res_merge$significance <- rowSums(cbind(res_merge$pFDR0<0.05, res_merge$pFDR1<0.05,
                                        res_merge$pFDR2<0.05, res_merge$pFDR4<0.05), na.rm = T)
addmargins(table(res_merge$Ann.,res_merge$pFDR0 < 0.05),1)
addmargins(table(res_merge$Ann.,res_merge$pFDR1 < 0.05),1)
addmargins(table(res_merge$Ann.,res_merge$pFDR2 < 0.05),1)
addmargins(table(res_merge$Ann.,res_merge$pFDR4 < 0.05),1)
addmargins(table(res_merge$Ann.,res_merge$significance),1)

write.csv(tab.res, file = "Result/Table 2/individual analysis_continuous/table2_ssf_lBMI.csv", row.names = F)
write.csv(res_merge, file = "Result/Table 2/individual analysis_continuous/table2_ssf_complete_lBMI.csv", row.names = F)

res_merge = res_merge %>% arrange(Ann.) %>% mutate(ID_new = 1:nrow(res_merge))
#library(ggplot2)
#library(ggrepel)

####--------------------- plot the overall beta plots for 4 visits ------------------- ####
res1 <- read.csv("Result/Table 2/individual analysis_continuous/bwz/table2_bwz_long_G1.csv")
res1 <- data.frame(res1, strata = "Male")
res2 <- read.csv("Result/Table 2/individual analysis_continuous/bwz/table2_bwz_long_G2.csv")
res2 <- data.frame(res2, strata = "Female")
res <- rbind(res1, res2)
res$Annotation <- factor(res$Annotation, levels = levels(forcats::fct_reorder(res$Annotation, as.numeric(res$Lipid.class))))
anno.tab <- unique(res[,c("Annotation","Lipid.class")])
anno.tab <- table(anno.tab$Annotation, anno.tab$Lipid.class)
anno.min <- (1:nrow(anno.tab))-0.5
anno.max <- (1:nrow(anno.tab))+0.5
min.list <- apply(anno.tab, 2, function(x) min(anno.min[x > 0]))
max.list <- apply(anno.tab, 2, function(x) max(anno.max[x > 0]))
y.anno <- data.frame(fill = as.factor(names(min.list)), ymin = min.list, ymax = max.list)
mycolors4 = brewer.pal(name="Dark2", n = 4)
mycolors12 = c(brewer.pal(name="Dark2", n = 6), brewer.pal(name="Paired", n = 6))
anno_break_sig <- unique(res$Annotation[res$pFDR < 0.05])
res$significance <- factor(res$pFDR < 0.05)
p <- ggplot(res, aes(x=beta, y = Annotation, color = factor(Visit)))+
  geom_pointrange(aes(xmin = beta - se, xmax = beta + se, shape = significance,alpha = as.numeric(significance)-0.7), position = ggstance::position_dodgev(height = 0.3))+
  scale_color_manual(name = "Visit", 
                     limits = c("0", "1", "2", "4"),
                     labels = c("0", "1", "2", "4"),
                     values = mycolors4)+
  scale_shape_manual(limits = c(F,T),
                     values = c(1,16),
                     guide = F)+
  scale_alpha_continuous(guide = F)+
  scale_y_discrete(name = "Annotation",
                   breaks = anno_break_sig,
                   labels = anno_break_sig,
                   guide = guide_axis(check.overlap = T))+
  theme(axis.text.y = element_text(size = 8))+
  theme_bw()+
  geom_vline(xintercept = 0) +
  geom_rect(y.anno, mapping = aes(ymin = ymin, ymax = ymax, fill=as.factor(fill)), 
            xmin = -Inf, xmax=Inf, alpha = 0.15, inherit.aes = F)+
  scale_fill_manual(name="Lipid Class",
                    limits = y.anno$fill,
                    labels = y.anno$fill,
                    values= mycolors12)+
  facet_wrap(~strata)

print(p)

####------------------------------ Volcano plots ----------------------------------
tit <- c("Birth weight","Birth-weight Z-score","Head Circumference",
         "Abdominal circumference","Length","Sum of Skinfolds")
ctt = 1
res_fin <- NULL
for(k in c("bw","bwz","head_circ","ab_circ","length","ssf")){
  res1 <- read.csv(paste0("Result/Table 2/individual analysis_continuous/",k,"/table2_",k,"_long_G1.csv"))
  res1$strata <-"Male"
  res2 <- read.csv(paste0("Result/Table 2/individual analysis_continuous/",k,"/table2_",k,"_long_G2.csv"))
  res2$strata <-"Female"
  res = rbind(res1,res2)
  res$label = res$Annotation
  res$significance <- factor(res$pFDR < 0.05,levels = c(FALSE,TRUE))
  res$label[res$significance==F] = NA
  res$tit <- tit[ctt]
  res_fin = rbind(res_fin,res)
  assign(paste0("p",ctt),
         ggplot(res,aes(x = beta, y = -log10(pFDR)))+
           geom_point(aes(col = as.factor(Lipid.class),shape = as.factor(Visit),alpha = as.numeric(significance)-0.6), size = 3)+
           geom_text_repel(aes(label = label),size = 3, max.overlaps = 10)+
           theme_bw()+scale_alpha_continuous(guide = F)+
           scale_color_manual(name="Lipid Class",
                              values= mycolors12)+
           scale_shape_discrete(name = "Visit") + 
           facet_wrap(~strata, scales = "free_x"))
  ctt = ctt + 1
}
ggplot(res_fin,aes(x = beta, y = -log10(pFDR)))+
  geom_point(aes(col = as.factor(Lipid.class), shape = as.factor(Visit), alpha = as.numeric(significance) - 0.6), size = 2)+
  geom_text_repel(aes(label = label), size =3,max.overlaps = 15)+
  theme_bw() + scale_alpha_continuous(guide = F)+
  scale_color_manual(name="Lipid Class",
                     values= mycolors12)+
  scale_shape_discrete(name = "Visit")+
  facet_wrap(strata~tit, scales = "free",nrow = 2)
for(i in 1:6){
  png(file = paste0("Figures/draft/visit-specific analysis_stratefied/Gender/",tit[i],"_volcano.png"),width = 1311,height=791)
  print(get(paste0("p",i)))
  dev.off()
}


#### Junk for plotting visit-specific plots
for(i in c(0,1,2,4)){
  p <- ggplot(res_merge, aes(x=ID_new, y= -log10(get(paste0("pFDR",i))), group=factor(Ann.),label = met_ID)) +
    geom_point(aes(col = Ann.,size = significance)) +
    scale_color_manual(values = mycolors)+
    scale_size(range = c(1.5,4), limits = c(0,4), guide = "none")+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    scale_y_continuous("-log10(FDR adjusted p-value)", 
                       limits=c(0, -log10(min(res$pFDR))))+
    scale_x_continuous("Metabolite ID", limits = c(1,328))+
    geom_hline(aes(yintercept = -log10(0.05)))+
    ggtitle(paste("Visit",i))+
    geom_text_repel(data = filter(res_merge, get(paste0("pFDR",i)) < 0.05), aes(label = met_ID), 
                    size = 3, box.padding = 0.5, max.overlaps = 10)
   pdf(file = paste0("Figures/individual analysis_continuous/ssf/pdf/ssf",i,"_lBMI.pdf"), wi = 8, he = 5)
     print(p)
   dev.off()
   png(file = paste0("Figures/individual analysis_continuous/ssf/ssf",i,"_lBMI.png"),width = 600,height=375)
   print(p)
   dev.off()
}

 for(i in c(0,1,2,4)){
    if(nrow(res_merge %>% filter(get(paste0("pFDR",i)) < 0.05)) > 0){
      p <- res_merge %>% filter(get(paste0("pFDR",i)) < 0.05) %>% ggplot(aes(x=get(paste0("beta",i)),y=ID_new,group = factor(Ann.))) + 
        geom_point(aes(col=factor(Ann.),size = significance)) + scale_color_manual(values = mycolors)+
        scale_size(range = c(1.5,4), limits = c(0,4), guide = "none")+
        guides(color = guide_legend(override.aes = list(size = 5)))+geom_vline(xintercept = 0) + 
        scale_y_discrete("Metabolite ID",breaks=NULL) + ggtitle(paste0("Visit",i)) + 
        scale_x_continuous(expression(beta),limits = c(min(res$beta),max(res$beta)))+
        geom_text_repel(data = filter(res_merge, get(paste0("pFDR",i)) < 0.05), aes(label = met_ID), 
                        size = 3, box.padding = 0.5)+
        labs(col = "Annotation")
      pdf(file = paste0("Figures/individual analysis_continuous/ssf/pdf/ssf",i,"_beta_lBMI.pdf"), wi = 8, he = 5)
      print(p)
      dev.off()
      png(file = paste0("Figures/individual analysis_continuous/ssf/ssf",i,"_beta_lBMI.png"),width = 600,height=375)
      print(p)
      dev.off()
  }else{next}
}

####---------- plot the betas for all 4 visits ---------####
res1 <- read.csv("Result/Table 2/individual analysis_continuous/table2_bw_long_G1.csv")
res2 <- read.csv("Result/Table 2/individual analysis_continuous/table2_bw_long_G2.csv")
res1 <- cbind(res1,gender = "Male")
res2 <- cbind(res2,gender = "Female")
if(!is.null(res1$convergence)){}
mycolors4 = brewer.pal(name="Dark2", n = 4)
names(res)
anno_lim <- c("1_SM 17:0 ","Acylcarnitine C18:2 ","CE (22:6) ","Ceramide (d42:1) ",
              "DG (38:6) ","Gal-Gal-Cer(d18:1/16:0) or Lactosylceramide(d18:1/16:0) ",
              "GlcCer (d42:2) ","Lactosylceramide (d18:1/24:1(15Z)) ", "LPC (p-18:0) or LPC (o-18:1) ",
              "PC (p-44:5) or PC (o-44:6) ", "PE (p-40:5) or PE (o-40:6) ", "SM (d44:2) ",
              "TG (60:2) ")
anno_labs <- c("1_SM", "Acylcarnitine", "CE", "Ceramide", "DG", "Gal", "GlcCer",
               "Lactosylceramide", "LPC","PC","PE","SM","TG")
anno_break_sig1 <- unique(res1$Annotation[res1$pFDR < 0.05])
anno_break_sig2 <- unique(res2$Annotation[res2$pFDR < 0.05])
res1$significance <- factor(res1$pFDR < 0.05)
res2$significance <- factor(res2$pFDR < 0.05)

p <- ggplot(res1, aes(x=beta, y = Annotation, color = factor(Visit)))+
  geom_pointrange(aes(xmin = beta - se, xmax = beta + se, shape = significance,alpha = as.numeric(significance)-0.7), position = ggstance::position_dodgev(height = 0.3))+
  scale_color_manual(name = "Visit", 
                     limits = c("0", "1", "2", "4"),
                     labels = c("0", "1", "2", "4"),
                     values = mycolors4)+
  scale_shape_manual(limits = c(F,T),
                     values = c(1,16),
                     guide = F)+
  scale_alpha_continuous(guide = F)+
  scale_y_discrete(name = "Annotation",
                   breaks = anno_break_sig1,
                   labels = anno_break_sig1,
                   guide = guide_axis(check.overlap = T))+
  theme(axis.text.y = element_text(size = 8))+
  geom_vline(xintercept = 0)+ggtitle("Male")

print(p)
