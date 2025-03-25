param <- commandArgs(trailingOnly=T)

POP  = eval(paste(text=param[1])) 
SNPset = eval(paste(text=param[2]))

data=read.table(paste0("TGP_PLINK.",SNPset,"_pca.no",POP,".eigenvec"),h=F)
rownames(data) = data$V2
INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise"
info=read.table(paste0(INFO_DIR, "/sorted_info.txt"),h=T,sep="\t")
mypop   = which(info$POP==POP)
mydata  = data[ mypop,]
dataTGP = data[-mypop,]
info    = info[-mypop,]
#compute center of each reference population
tokeep=7
center_afr = apply(subset(dataTGP[,3:tokeep],info$myREG=="Africa"),2,mean)
center_amr = apply(subset(dataTGP[,3:tokeep],info$myREG=="American"),2,mean)
center_eas = apply(subset(dataTGP[,3:tokeep],info$myREG=="East Asia"),2,mean)
center_eur = apply(subset(dataTGP[,3:tokeep],info$myREG=="European"),2,mean)
center_mea = apply(subset(dataTGP[,3:tokeep],info$myREG=="Middle-East"),2,mean)
center_sas = apply(subset(dataTGP[,3:tokeep],info$myREG=="South Asia"),2,mean)
out=NULL
for (i in 1:nrow(mydata)){
    thisiid = c(
    sum((mydata[i,3:tokeep]-center_afr)**2),
    sum((mydata[i,3:tokeep]-center_amr)**2),
    sum((mydata[i,3:tokeep]-center_eas)**2),
    sum((mydata[i,3:tokeep]-center_eur)**2),
    sum((mydata[i,3:tokeep]-center_mea)**2),
    sum((mydata[i,3:tokeep]-center_sas)**2))
    thisiid = as.integer(thisiid == min(thisiid))
    out=rbind(out,thisiid)
}
colnames(out) = c("afr","amr","eas","eur","mea","sas")
rownames(out)= mydata[,2]
write.table(out,file=paste0("pca_results/",SNPset,".no",POP,".txt"),col.names=T,row.names=F,quote=F,sep="\t")

#
library(randomForest)
ancestry=as.character(info$myREG)
dataRF = cbind(ancestry,dataTGP[,3:tokeep])
model <- randomForest(as.factor(ancestry) ~ ., data=dataRF)
pred  <- predict(model, mydata[,3:tokeep], type="prob")
write.table(pred,file=paste0("pcaRF_results/",SNPset,".no",POP,".txt"),col.names=T,row.names=F,quote=F,sep="\t")



