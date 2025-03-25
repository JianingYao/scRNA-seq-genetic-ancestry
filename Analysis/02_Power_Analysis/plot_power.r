param <- commandArgs(trailingOnly=T)

STUDY = eval(paste(text=param[1]))
SNPset = eval(paste(text=param[2])) 

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis"
setwd(paste0(POWER_DIR, "/Results/", STUDY))
INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")[,3:4]
if (length(which(duplicated(info))) > 0) {
  info = info[-which(duplicated(info)),]
}

pca.popres  = NULL
pcaRF.popres = NULL
adm.popres = NULL
for (i in (1:nrow(info))){
    pop = info$POP[i]
    reg = info$REG[i]
    #
    pca = read.table(paste0(POWER_DIR, "/Results/", STUDY, "/pca_results/", SNPset, ".no",pop,".txt"),h=T)
    colnames(pca) = c("afr", "amr", "eas", "eur", "mid", "sas")
    pca = pca[,as.character(reg)]
    pca.popres = rbind(pca.popres,c(sum(pca),length(pca)))
    #
    pcaRF = read.table(paste0(POWER_DIR, "/Results/", STUDY, "/pcaRF_results/", SNPset, ".no",pop,".txt"),h=F, skip=1)
    colnames(pcaRF) = c("afr", "amr", "eas", "eur", "mid", "sas")
    mymax = apply(pcaRF,1,max)
    pcaRF = pcaRF[,as.character(reg)]
    pcaRF.popres = rbind(pcaRF.popres,c(sum(pcaRF==mymax),length(pcaRF)))
    #
    adm = read.table(paste0(POWER_DIR, "/Results/", STUDY, "/admixture_results/", SNPset, ".no",pop,".txt"),h=F)
    colnames(adm) = unique(info$REG[-i])
    mymax = apply(adm,1,max)
    adm = adm[,as.character(reg)]
    adm.popres = rbind(adm.popres,c(sum(adm==mymax),length(adm)))
}
# overall error rate
pca.popres = cbind(info,pca.popres)
pcaRF.popres = cbind(info,pcaRF.popres)
adm.popres = cbind(info,adm.popres)
pca_accuracy = round(100*(1 - sum(pca.popres[,3])/sum(pca.popres[,4])),2)
rf_accuracy = round(100*(1 - sum(pcaRF.popres[,3])/sum(pcaRF.popres[,4])),2)
adm_accuracy = round(100*(1 - sum(adm.popres[,3])/sum(adm.popres[,4])),2)
accuracy = cbind(pca_accuracy, rf_accuracy, adm_accuracy)
colnames(accuracy) = c("PCA(%)", "Random Forest(%)", "ADMIXTURE(%)")
write.table(accuracy,file=paste0(POWER_DIR, "/Results/", STUDY, "/", STUDY, "_error.txt"),col.names=T,row.names=F,quote=F,sep="\t")

# continent-wise error rate
pca.regres = NULL
adm.regres = NULL
pcaRF.regres = NULL
for (reg in unique(info$REG)){
    pca.regres = rbind(pca.regres, colSums(subset(pca.popres[,3:4], pca.popres[,2]==reg)))
    adm.regres = rbind(adm.regres, colSums(subset(adm.popres[,3:4], adm.popres[,2]==reg)))
    pcaRF.regres = rbind(pcaRF.regres, colSums(subset(pcaRF.popres[,3:4], pcaRF.popres[,2]==reg)))
}

pca.regres = cbind(pca.regres, round(100*(1 - pca.regres[,1]/pca.regres[,2]),2))
pca.regres = as.data.frame(pca.regres)
rownames(pca.regres) = unique(info$REG)
colnames(pca.regres) = c('pred', 'reference', 'error(%)') 
write.table(pca.regres,file=paste0(POWER_DIR, "/Results/", STUDY, "/", STUDY, "_pca_continent_error.txt"),col.names=T,row.names=T,quote=F,sep="\t")

pcaRF.regres = cbind(pcaRF.regres, round(100*(1 - pcaRF.regres[,1]/pcaRF.regres[,2]),2))
pcaRF.regres = as.data.frame(pcaRF.regres)
rownames(pcaRF.regres) = unique(info$REG)
colnames(pcaRF.regres) = c('pred', 'reference', 'error(%)') 
write.table(pcaRF.regres,file=paste0(POWER_DIR, "/Results/", STUDY, "/", STUDY, "_pcaRF_continent_error.txt"),col.names=T,row.names=T,quote=F,sep="\t")

adm.regres = cbind(adm.regres, round(100*(1 - adm.regres[,1]/adm.regres[,2]),2))
adm.regres = as.data.frame(adm.regres)
rownames(adm.regres) = unique(info$REG)
colnames(adm.regres) = c('pred', 'reference', 'error(%)') 
write.table(adm.regres,file=paste0(POWER_DIR, "/Results/", STUDY, "/", STUDY, "_admixture_continent_error.txt"),col.names=T,row.names=T,quote=F,sep="\t")


