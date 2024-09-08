library(sva)
library(pvca)
library(umap)

library(reshape2)
lipid <- read.csv("Data/lipid_imputed_wosmoke.csv") 
met_idx <- 4:423
DFmet <- t(lipid[,met_idx])
mod = model.matrix(~GA_blood,data=lipid)
batch = lipid$batch

loglipid <- lipid
loglipid[,met_idx] = log(lipid[,met_idx])
loglipid.data = loglipid[,met_idx]
loglipid.labels = loglipid$batch
lipid.umap = umap(loglipid.data)
head(lipid.umap$layout,3)
plot.lipid(lipid.umap, loglipid.labels)

combat_DFmet = ComBat(dat= log(DFmet), batch = batch, mod = mod, par.prior = T, prior.plots = T)
lipid_new <- lipid
lipid_new[,met_idx] <- t(combat_DFmet)
write.csv(lipid_new, "Data/lipid_imputed_wosmoke_combat.csv", row.names = F)

lipid <- read.csv("../Data/lipid_imputed_wosmoke_combat.csv")
lipid$fasting <- factor(lipid$diet > 8)
# calculate percentage of fasting at each visit
print.table(round(table(lipid$Visit,lipid$fasting) / addmargins(table(lipid$Visit,lipid$fasting),2)[,3]*100,1),
       row.names = F)

lipid.combat.data <- lipid[,met_idx]
lipid.combat.labels <- lipid$batch
lipid.combat.umap <- umap(lipid.combat.data)
plot.lipid(lipid.combat.umap, lipid.combat.labels)

lipid.combat.labels <- lipid$fasting
lipid.combat.umap <- umap(lipid.combat.data[lipid$Visit==4,])
par(mfrow=c(2,2))
plot.lipid(lipid.combat.umap, lipid.combat.labels[lipid$Visit==4])
lipid
lipid_wt <- lipid %>% select(ID,wt)
lipid_wt <- aggregate(wt~ID, data = lipid_wt, head, 1)

lipid_comp_wt <- lipid %>% select(ID,Visit,wt) %>%filter(Visit %in% c(2,4))
lipid_comp_wt <- aggregate(wt~ID, data = lipid_comp_wt, head, 1)
met_known <- 4:331
met_all <- 4:423

data_time <- lipid[,c(1,3,2)]
names(data_time) <- c("ID", "Visit", "Time")
data_time <- dcast(data_time, ID ~ Visit)
names(data_time)[2:5] <- paste0("T",1:4)

data_diet <- lipid[,c(1,3,463)]
names(data_diet) <- c("ID", "Visit", "Diet")
data_diet[3] <- ifelse(data_diet[3] <= 8, 0, 1)
data_diet <- dcast(data_diet, ID ~ Visit)
names(data_diet)[2:5] <- paste0("TCOV",1:4)

data_cov <- inner_join(data_time, data_diet, by="ID")

for(i in met_all){
  datai_met <- lipid[,c(1,3,i)]
  names(datai_met) <- c("ID", "Visit", "Met")
  datai_met <- dcast(datai_met, ID ~ Visit)
  names(datai_met)[2:5] <- paste0("V", 1:4)
  
  datai <- merge(data_cov, datai_met, by = "ID")
  datai <- merge(datai,lipid_wt, by = "ID", all.x = T)
  write.csv(datai, paste0("Data/sas_data/data",i-3,".csv"), row.names = F)
}

for(i in met_all){
  datai_met <- lipid[,c(1,3,i)]
  names(datai_met) <- c("ID", "Visit", "Met")
  datai_met <- dcast(datai_met, ID ~ Visit)
  names(datai_met)[2:5] <- paste0("V", 1:4)
  
  datai_comp <- merge(data_cov, datai_met, by = "ID")
  datai_comp <- inner_join(lipid_comp_wt, datai_comp, by = "ID")
  write.csv(datai_comp, paste0("Data/sas_data_comp/data",i-3,".csv"), row.names = F)
}

plot.lipid <- function(x,labels, main="A UMAP visualization of the lipid(log-tranformed) dataset",
                       colors=c("#ff7f00", "#e377c2", "#17becf"),
                       pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                       cex.main=1, cex.legend=0.85) {
    layout = x
    if (is(x, "umap")) {
      layout = x$layout
    }

    xylim = range(layout)
    xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
    if (!add) {
      par(mar=c(0.2,0.7,1.2,0.7), ps=10)
      plot(xylim, xylim, type="n", axes=F, frame=F)
      rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
    }
    points(layout[,1], layout[,2], col=colors[as.integer(labels)],
           cex=cex, pch=pch)
    mtext(side=3, main, cex=cex.main)

    labels.u = unique(labels)
    legend.pos = "topleft"
    legend.text = as.character(labels.u)
    if (add) {
      legend.pos = "bottomleft"
      legend.text = paste(as.character(labels.u), legend.suffix)
    }

    legend(legend.pos, legend=legend.text, inset=0.03,
           col=colors[as.integer(labels.u)],
           bty="n", pch=pch, cex=cex.legend)
}
