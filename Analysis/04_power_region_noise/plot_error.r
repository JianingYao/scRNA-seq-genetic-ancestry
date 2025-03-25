param <- commandArgs(trailingOnly=T)
ERROR = as.character(eval(paste(text=param[1]))) 

library(ggplot2)
library(data.table)

SNPset = c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia", "eQTLAutoimmune")

# admixture
adm_list <- lapply(SNPset, function(SNP) {
    read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise/Results_", ERROR, "/", SNP, "/final_outs/", SNP, "_admixture_anc_accuracy.txt"), header = TRUE, col.names = c(SNP))
})
adm_error <- 100*(1 - do.call(cbind, adm_list))
write.table(adm_error, file=paste0("./Results_", ERROR, "/admixture_error.txt"), col.names=T,row.names=T,quote=F,sep="\t")
adm_error$Continent <- rownames(adm_error)
adm_long <- melt(setDT(adm_error), id.vars = "Continent", variable.name = "SNP", value.name = "Error")
adm_plot <- ggplot(adm_long, aes(x = Continent, y = Error, fill = SNP)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Continent", y = "Error (%)", title = "Admixture error rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave(paste0("./Results_", ERROR, "/admixture_error.png"), adm_plot, width = 10, height = 6, dpi = 300)

# pca
pca_list <- lapply(SNPset, function(SNP) {
    read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise/Results_", ERROR, "/", SNP, "/final_outs/", SNP, "_pca_anc_accuracy.txt"), header = TRUE, col.names = c(SNP))
})
pca_error <- 100*(1 - do.call(cbind, pca_list))
write.table(pca_error, file=paste0("./Results_", ERROR, "/pca_error.txt"), col.names=T,row.names=T,quote=F,sep="\t")
pca_error$Continent <- rownames(pca_error)
pca_long <- melt(setDT(pca_error), id.vars = "Continent", variable.name = "SNP", value.name = "Error")
pca_plot <- ggplot(pca_long, aes(x = Continent, y = Error, fill = SNP)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Continent", y = "Error (%)", title = "PCA error rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave(paste0("./Results_", ERROR, "/pca_error.png"), pca_plot, width = 10, height = 6, dpi = 300)

# rf
rf_list <- lapply(SNPset, function(SNP) {
    read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise/Results_", ERROR, "/", SNP, "/final_outs/", SNP, "_RF_anc_accuracy.txt"), header = TRUE, col.names = c(SNP))
})
rf_error <- 100*(1 - do.call(cbind, rf_list))
write.table(rf_error, file=paste0("./Results_", ERROR, "/rf_error.txt"), col.names=T,row.names=T,quote=F,sep="\t")
rf_error$Continent <- rownames(rf_error)
rf_long <- melt(setDT(rf_error), id.vars = "Continent", variable.name = "SNP", value.name = "Error")
rf_plot <- ggplot(rf_long, aes(x = Continent, y = Error, fill = SNP)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Continent", y = "Error (%)", title = "RF error rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave(paste0("./Results_", ERROR, "/rf_error.png"), rf_plot, width = 10, height = 6, dpi = 300)
