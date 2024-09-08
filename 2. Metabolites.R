##############################################################################
## Single metabolite analysis, adjusting batch effect via meta-analysis
##############################################################################
rm(list = ls())
library(nlme)
library(rmeta)
#setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
anno <- read.csv("../Data/Lipidomics_annotation_file.csv",header=T, sep=",")
anno$met_ID <- chartr(old = ".", new = "_", anno$lipid_identifier)

for (visit in c(0, 1, 2, 4)) {    # By visit
 
  ###1.Read and prepare data---------------------------------------------------------------
  DF <- read.csv(paste("../tmp/lipid_Std", visit, ".csv", sep=""), header=T, sep=",")
  DF <- data.frame(ID = DF$ID, pairID = DF$PairID, v = (1/DF$wt)^2, GA_blood = DF$GA_blood,
                   age = DF$age, race = factor(DF$race), education = factor(DF$Education, ordered = T), 
                   preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = DF$Gender,
                   anthro = DF$Circ, alcohol = DF$Alcohol_bin, batch = factor(DF$batch), post = DF$measuredate,
                   DF[,2:421])
  met_index <- c(15:434)
  
  ###2.Estimate batch-specific effects-----------------------------------------------------
  #Define LME function: random intercept - pairID, weighted
  lme_met <- function(x) {
    nDF2 <- nDF[,c(1:14)]
    nDF2$met <- x
    nDF2 <- na.omit(nDF2)
    #model = lme(anthro ~ met + age + as.factor(race) + preBMI + parity + female + GADel + glucose, 
    #            data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    model = lme(anthro ~ met + age + race + education + preBMI + GA_blood + Gender + nulli + alcohol + post, 
                data = nDF2, random = ~1 | pairID, weights = varFixed(~v))
    est   <- c(summary(model)$tTable[2,])
    return(est)
  }
  
  #Estimate by batch
  for (k in 1:3) {
    nDF <- DF[DF$batch == k,]
    res <- as.data.frame(t(sapply(nDF[,met_index], lme_met)))[,c(1,2,5)]
    colnames(res) <- c("beta","se","pval")
    res$met_ID <- substring(row.names(res), 4)
    res <- merge(anno, res, by = "met_ID", all.x=F, all.y=T)
    write.csv(res, paste("tmp/LME_met", visit, "_batch", k, ".csv",sep=""), row.names=F)
  }
  ###3.Meta-analysis across batch: inverse variance weight--------------------------------
  
  res1 <- read.csv(paste("tmp/LME_met", visit, "_batch1.csv",sep=""), header=T, sep=",")
  res2 <- read.csv(paste("tmp/LME_met", visit, "_batch2.csv",sep=""), header=T, sep=",")
  res3 <- read.csv(paste("tmp/LME_met", visit, "_batch3.csv",sep=""), header=T, sep=",")
  
  
  meta <- data.frame(lipid_id = res1$lipid_identifier, annotation = res1$Annotation, lb  = rep(NA, nrow(res1)),
                          beta = rep(NA,nrow(res1)), ub = rep(NA,nrow(res1)), pval = rep(NA,nrow(res1)))
  for(i in 1:420){
    res = meta.summaries(c(res1$beta[i],res2$beta[i],res3$beta[i]), c(res1$se[i],res2$se[i],res3$se[i]), method = "random")
    meta[i,3:5] = summary(res)$summci
    meta[i,6] = res$test[2]
  }
  meta$pBonf = p.adjust(meta$pval,"bonferroni")
  meta$pFDR  = p.adjust(meta$pval,"fdr")
  
  meta <- meta[order(meta$pval),]
  write.csv(meta, paste("tmp/LME_met", visit, "_meta.csv",sep=""), row.names=F)
  
  ###4.QQ plot-----------------------------------------------------------------------------
  #png(paste("QQplot_visit", visit, ".png", sep=""), width=2400,height=2400,res=500)
  #pval = meta$pval
  #exp = -log10(rank(pval)/(length(pval)+1))
  #obs = -log10(pval)
  #plot(exp,obs,xlab="Expected -log10(P)",ylab="Observed -log10(P)",
  #     xlim=c(0,6),ylim=c(0,6))
  #abline(0,1)
  #title(paste("Lipidomics metabolites: meta p-value, visit ", visit))
  #dev.off()

}

#Merge results from all 4 visits
meta0 <- read.csv("tmp/LME_met0_meta.csv", header=T, sep=",")
meta1 <- read.csv("tmp/LME_met1_meta.csv", header=T, sep=",")
meta2 <- read.csv("tmp/LME_met2_meta.csv", header=T, sep=",")
meta4 <- read.csv("tmp/LME_met4_meta.csv", header=T, sep=",")
for(visit in c(0,1,2,4)){
  assign(paste0("meta",visit), data.frame(get(paste0("meta",visit)), Visit = visit))
}

meta = rbind(meta0,meta1,meta2,meta4)


tab.res = meta %>% filter(pFDR < 0.05) %>%
  arrange(annotation, lipid_id) %>%
  select(lipid_id,annotation,Visit,beta, lb, ub)

write.csv(tab.res, "Result/_LME_met_visits_meta_circ.csv", row.names=F)

# meta.all <- Reduce(function(x,y) merge(x = x, y = y, by = "lipid_identifier"), 
#                    list(meta0, meta1, meta2, meta4))
# write.csv(meta.all, "Result/_LME_met_visits_meta_circ.csv", row.names=F)

# Obtain Table 2 by only keeping the significant results.
bw <- read.csv("Result/_LME_met_visits_meta_bw1000.csv")
bwz <- read.csv("Result/_LME_met_visits_meta_bwz.csv")
length <- read.csv("Result/_LME_met_visits_meta_length.csv")
circ <- read.csv("Result/_LME_met_visits_meta_circ.csv")
#circ$Ann. <- factor(substr(circ$Annotation,1,3))
#table(circ$Ann., circ$Visit, circ$beta < 0)[,,1]
ssf <- read.csv("Result/_LME_met_visits_meta_ssf.csv")


all <- Reduce(function(x,y) merge(x = x, y = y, by = c("lipid_id","annotation","Visit"), all = T), 
                                  list(bw,bwz,length,circ,ssf))
all$beta.bw <- paste0(round(all[,4],2), " (", round(all[,5],2), ", ",round(all[,6],2),")")
all$beta.bw[is.na(all[,4])] <- NA
all$beta.bwz <- paste0(round(all[,7],2), " (", round(all[,8],2), ", ",round(all[,9],2),")")
all$beta.bwz[is.na(all[,7])] <- NA
all$beta.length <- paste0(round(all[,10],2), " (", round(all[,11],2), ", ",round(all[,12],2),")")
all$beta.length[is.na(all[,10])] <- NA
all$beta.circ <- paste0(round(all[,13],2), " (", round(all[,14],2), ", ",round(all[,15],2),")")
all$beta.circ[is.na(all[,13])] <- NA
all$beta.ssf <- paste0(round(all[,16],2), " (", round(all[,17],2), ", ",round(all[,18],2),")")
all$beta.ssf[is.na(all[,16])] <- NA

library(ggplot2)
library("RColorBrewer")
library(gtools)
library(pals)

all.plot <- all[,c(2,1,3,4,7,10,13,16)]
names(all.plot)[c(4:8)] <- c("bw","bwz","length","circ","ssf")
all.plot$Ann. <- substring(all.plot$annotation,1,3)

all.plot %>%
  arrange(Ann., lipid_id)%>%
  ggplot(aes(x = circ, y = lipid_id))+
  geom_point(aes(colour = Ann., shape=Ann.), size = 3)+
  scale_shape_manual(values = rep(c(8, 16, 17),each = 5))+
  scale_colour_manual(values = rep(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF"),3), na.value="white")+
  geom_vline(xintercept = 0)+
  labs(x = expression(paste(beta," for circumferences (mm)")),y = "Lipid")+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),)+
  facet_grid(cols = vars(Visit))+
  theme(panel.background = element_rect(fill="white",colour = "black"))

all <- all[c(2,1,3,19:23)]
all$annotation <- factor(all$annotation,levels = sort(levels(all$annotation)))
all <- all %>% arrange(annotation,lipid_id,Visit)
write.csv(all, "Result/Table2.csv", row.names = F)
