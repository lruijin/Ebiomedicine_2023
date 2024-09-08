load("Visit-0-networkConstruction-man.RData")

pdf("~/Collaboration_fetal/lipidomics/version 6.0/Figure4_a.pdf",
    width = 3.02, height = 2.24)
par(cex=0.4)
plot(METree, main = "Clustering of modules at Visit 0", xlab = "", sub = "")
dev.off()
load("Visit-1-networkConstruction-man.RData")
pdf("~/Collaboration_fetal/lipidomics/version 6.0/Figure4_b.pdf",
    width = 3.02, height = 2.24)
par(cex=0.4)
plot(METree, main = "Clustering of modules at Visit 1", xlab = "", sub = "")
dev.off()
load("Visit-2-networkConstruction-man.RData")
pdf("~/Collaboration_fetal/lipidomics/version 6.0/Figure4_c.pdf",
    width = 3.02, height = 2.24)
par(cex=0.4)
plot(METree, main = "Clustering of modules at Visit 2", xlab = "", sub = "")
dev.off()
load("Visit-4-networkConstruction-man.RData")
pdf("~/Collaboration_fetal/lipidomics/version 6.0/Figure4_d.pdf",
    width = 3.02, height = 2.24)
par(cex=0.4)
plot(METree, main = "Clustering of modules at Visit 4", xlab = "", sub = "")
dev.off()
