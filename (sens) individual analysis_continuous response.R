setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
rm(list= ls())
lipid <- read.csv("../Data/lipid_imputed_wosmoke_combat.csv", header = T)
# anno <- read.csv("../Data/Lipidomics_annotation_file.csv",header=T, sep=",")
# anno_class <- read.csv("../Data/_Lipid class.csv", header = T, sep=",")
# anno <- merge(anno_class, anno, by="lipid_identifier", all.x = T)
library(ggplot2)

# check if the two files match to each other
# sum(as.character(anno$Annotation.x) != as.character(anno$Annotation.y))
# make a column for complete annotation to add m/z and rt to the replicated annotations
# Annotation <- as.character(anno$Annotation.x)
# Annotation_comp <- paste0(anno$Annotation.x," (m/z: ", anno$Batch.m.z, ", rt: ",anno$Batch.RT,")")
# Annotation[duplicated(Annotation)|duplicated(Annotation, fromLast= T)] <- 
#     Annotation_comp[duplicated(Annotation) | duplicated(Annotation,fromLast=T)]
# anno$Annotation_comp <- Annotation
# anno$Annotation <- anno$Annotation.x
# anno <- anno %>% select(-c(Annotation.x,Annotation.y)) %>%
#     mutate(met_ID = chartr(old = ".", new = "_", anno$lipid_identifier))

# write.csv(anno, "Data/lipidomics_annotation_complete.csv",row.names = F)
anno <- read.csv("../Data/lipidomics_annotation_complete.csv")
library(lme4)
library(lmerTest)
library(dplyr)
met_idx = 4:423
anthros <- c("bw","bwz","Length","ab_circ","head_circ","SSF")
for (visit in c(0, 1, 2, 4)) {    # By visit
    for(j in 1:length(anthros)){
        ###1.Read and prepare data---------------------------------------------------------------
        DF <- lipid[lipid$Visit == visit,]
        if(j == 1){
            DF <- data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                             age = DF$age, race = factor(DF$race==3), 
                             education = factor(DF$Education > 3, ordered = T), 
                             preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), 
                             Gender = factor(DF$Gender), 
                             GWG = DF[,paste0("GWGv",visit)],
                             anthro = DF$bw.x/1000, DF[,met_idx])
        }else{
            anthro_idx = which(names(DF) == anthros[j])
            DF <- data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                             age = DF$age, race = factor(DF$race==3), 
                             education = factor(DF$Education > 3, ordered = T), 
                             preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), 
                             Gender = factor(DF$Gender),
                             GWG = DF[,paste0("GWGv",visit)],
                             anthro = DF[,anthro_idx], DF[,met_idx])
        }
        
        # met_index_df <- c(15:434) # all lipids
        met_index_df <- c(13:340) # all known lipids 328 intotal
        
        ###2.Estimate batch-specific effects-----------------------------------------------------
        #Define LME function: random intercept - pairID, weighted
        lme_met <- function(x) {
            DF2 <- DF[,c(1:12)]
            DF2$met <- x
            DF2 <- na.omit(DF2)
            model = lmer(anthro ~ met + age + race + education + preBMI + GA_blood + Gender + GWG + nulli + (1 | pairID), 
                         data = DF2, weights = 1/v)
            #model = lme(anthro ~ met + age + race + education + preBMI + GA_blood + Gender + nulli, 
            #            data = DF2, random = ~1 | pairID, weights = varFixed(~v))
            est   <- c(summary(model)$coefficients[2,])
            return(est)
        }
        
        
        res <- as.data.frame(t(sapply(DF[,met_index_df], lme_met)))[,c(1,2,5)]
        colnames(res) <- c("beta","se","pval")
        res$met_ID <- substring(row.names(res), 4)
        res <- inner_join(res,anno, by = "met_ID")
        res$pBonf <- p.adjust(res$pval,"bonferroni")
        res$pFDR <- p.adjust(res$pval,"fdr")
        
        write.csv(res, paste("../tmp/known/", anthros[j], visit, ".csv",sep=""), row.names=F)
    }
}

## ------------------- Part 2: summary the result ------------------------------##
for(j in 1:length(anthros)){
    for(i in c(0,1,2,4)){
        assign(paste0("res",i), read.csv(paste0("../tmp/known/",anthros[j],i,".csv"), header=T, sep=","))
        assign(paste0("res",i), data.frame(get(paste0("res",i)), Visit = i))
    }
    
    res = rbind(res0,res1,res2,res4)
    
    
    tab.res = res %>% filter(pFDR < 0.05) %>%
        arrange(Lipid.class, met_ID) %>%
        select(Lipid.class,Annotation,Annotation_comp, met_ID,Visit,beta,pval,pFDR)
    
    res_merge <- inner_join(inner_join(inner_join(res0[,c("Lipid.class","Annotation", "Annotation_comp", "met_ID","beta","se","pval","pBonf","pFDR")],
                                                  res1[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID"),
                                       res2[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID"),
                            res4[,c("met_ID","beta","se","pval","pBonf","pFDR")], by="met_ID")
    
    names(res_merge)[c(5:24)] <- paste0(rep(c("beta","se","pval","pBonf","pFDR"),4),rep(c(0,1,2,4),each=5))
    res_merge$significance <- (res_merge$pFDR0<0.05) + (res_merge$pFDR1<0.05) + (res_merge$pFDR2<0.05) + (res_merge$pFDR4<0.05)
    # addmargins(table(res_merge$Ann.,res_merge$pFDR0 < 0.05),1)
    # addmargins(table(res_merge$Ann.,res_merge$pFDR1 < 0.05),1)
    # addmargins(table(res_merge$Ann.,res_merge$pFDR2 < 0.05),1)
    # addmargins(table(res_merge$Ann.,res_merge$pFDR4 < 0.05),1)
    # addmargins(table(res_merge$Ann.,res_merge$significance),1)
    write.csv(res,file = paste0("../Result/revision/individual analysis_continuous/table2_",anthros[j],"_long.csv"), row.names = F)
    write.csv(tab.res, file = paste0("../Result/revision/individual analysis_continuous/table2_",anthros[j],".csv"), row.names = F)
    write.csv(res_merge, file = paste0("../Result/revision/individual analysis_continuous/table2_",anthros[j],"_complete.csv"), row.names = F)
    
}

####--------------------- plot the overall beta plots for 4 visits ------------------- ####
library(RColorBrewer)
library(ggplot2)
res <- read.csv("../Result/revision/individual analysis_continuous/table2_SSF_long.csv")
res$Annotation <- factor(res$Annotation, levels = levels(forcats::fct_reorder(res$Annotation, as.numeric(factor(res$Lipid.class)))))

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
    geom_pointrange(aes(xmin = beta - se, xmax = beta + se, shape = significance,
                        alpha = as.numeric(significance)-0.7), 
                    position = ggstance::position_dodgev(height = 0.3))+
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
    geom_vline(xintercept = 0) +
    geom_rect(y.anno, mapping = aes(ymin = ymin, ymax = ymax, fill=as.factor(fill)), 
              xmin = -Inf, xmax=Inf, alpha = 0.15, inherit.aes = F)+
    scale_fill_manual(name="Lipid Class",
                      limits = y.anno$fill,
                      labels = y.anno$fill,
                      values= mycolors12)

ggsave(p, filename = "../Figures/revision/visit-specific analysis_mixed effect model/ssf.png",
       height = 13)
####------------------------------ Volcano plots ----------------------------------
library(ggrepel)
tit <- c("Birth weight","Birth-weight Z-score","Head Circumference",
         "Abdominal circumference","Length","Sum of Skinfolds")
ctt = 1
res_fin <- NULL
for(k in c("bw","bwz","head_circ","ab_circ","length","ssf")){
    res <- read.csv(paste0("../Result/revision/individual analysis_continuous/table2_",k,"_long.csv"))
    res$label = res$Annotation
    res$significance <- factor(res$pFDR < 0.05,levels = c(FALSE,TRUE))
    res$label[res$significance==F] = NA
    res$tit <- tit[ctt]
    res_fin = rbind(res_fin,res)
    assign(paste0("p",ctt),
           ggplot(res,aes(x = beta, y = -log10(pFDR)))+
               geom_point(aes(col = as.factor(Lipid.class),shape = as.factor(Visit),alpha = as.numeric(significance)-0.6), size = 3)+
               geom_text_repel(aes(label = label),size = 3, max.overlaps = 30)+
               theme_bw()+scale_alpha_continuous(guide = F)+
               scale_color_manual(name="Lipid Class",
                                  values= mycolors12)+
               scale_shape_discrete(name = "Visit"))
    ctt = ctt + 1
}
ggplot(res_fin,aes(x = beta, y = -log10(pFDR)))+
    geom_point(aes(col = as.factor(Lipid.class), shape = as.factor(Visit), alpha = as.numeric(significance) - 0.6), size = 2)+
    geom_text_repel(aes(label = label), size =3,max.overlaps = 15)+
    theme_bw() + scale_alpha_continuous(guide = F)+
    scale_color_manual(name="Lipid Class",
                       values= mycolors12)+
    scale_shape_discrete(name = "Visit")+
    facet_wrap(~tit, scales = "free_x")

ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3,labels = tit, common.legend = T,
          legend = "bottom", vjust = -0.1)+
    theme(plot.margin = margin(t=1.2,r = 0.4, b = 0.4, l =0.4, unit = "cm"))
ggsave(filename = "../Figures/revision/visit-specific analysis_mixed effect model/volcano.png")
#### Junk for plotting visit-specific plots
res_merge = res_merge %>% arrange(Lipid.class) %>% mutate(ID_new = 1:nrow(res_merge))
#library(ggplot2)
#library(ggrepel)
#library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
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
                        size = 3, box.padding = 0.5, max.overlaps = 50)
    #pdf(file = paste0("Figures/individual analysis_continuous/",anthros[j],"/pdf/", anthros[j] ,i,".pdf"), wi = 8, he = 5)
    #  print(p)
    #dev.off()
    png(file = paste0("Figures/individual analysis_continuous/",anthros[j],"/", anthros[j] ,i,".png"),width = 600,height=375)
    print(p)
    dev.off()
}

for(i in c(0,1,2,4)){
    if(nrow(res_merge %>% filter(get(paste0("pFDR",i)) < 0.05)) > 0){
        p <- res_merge %>% filter(get(paste0("pFDR",i)) < 0.05) %>% 
            ggplot(aes(x=get(paste0("beta",i)),y=ID_new,group = factor(Ann.))) + 
            geom_point(aes(col=factor(Ann.),size = significance)) + scale_color_manual(values = mycolors)+
            scale_size(range = c(1.5,4), limits = c(0,4), guide = "none")+
            guides(color = guide_legend(override.aes = list(size = 5)))+geom_vline(xintercept = 0) + 
            scale_y_discrete("Metabolite ID",breaks=NULL) + ggtitle(paste0("Visit",i)) + 
            scale_x_continuous(expression(beta),limits = c(min(res$beta),max(res$beta)))+
            geom_text_repel(data = filter(res_merge, get(paste0("pFDR",i)) < 0.05), aes(label = met_ID), 
                            size = 3, box.padding = 0.5,max.overlaps = 50)+
            labs(col = "Annotation")
        #pdf(file = paste0("Figures/individual analysis_continuous/",anthros[j],"/pdf/",anthros[j],"_beta",i,".pdf"), wi = 8, he = 5)
        #print(p)
        #dev.off()
        png(file = paste0("Figures/individual analysis_continuous/",anthros[j],"/",anthros[j],"_beta",i,".png"), width = 600, height = 375)
        print(p)
        dev.off()
    }else{next}
}