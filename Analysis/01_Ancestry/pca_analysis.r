param <- commandArgs(trailingOnly=T)
STUDY = eval(paste(text=param[1])) 
DIR_NAME = eval(paste(text=param[2]))

library(RColorBrewer)

DATA_DIR=paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/", DIR_NAME)
INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
data=read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".TGP_HGDP.pca.eigenvec"),h=F)
rownames(data) = data$V2
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
rownames(info) = info$IID
mydata  = data[as.character(read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".plink.fam"))[,2]),]
dataTGP = data[as.character(info$IID),]

topPCtokeep = 5  #let's keep the 5 top PCs
tokeep = topPCtokeep+2

#compute center of each reference population
center_afr = apply(subset(dataTGP[,3:tokeep],info$myREG=="Africa"),2,mean)
center_amr = apply(subset(dataTGP[,3:tokeep],info$myREG=="American"),2,mean)
center_eas = apply(subset(dataTGP[,3:tokeep],info$myREG=="East Asia"),2,mean)
center_eur = apply(subset(dataTGP[,3:tokeep],info$myREG=="European" | info$myREG=="Finnish"),2,mean)
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
if (!dir.exists(paste0(DATA_DIR,"/", STUDY,"/PCA_results"))) {
  dir.create(paste0(DATA_DIR,"/", STUDY,"/PCA_results"))
}
write.table(out,file=paste0(DATA_DIR,"/", STUDY,"/PCA_results/",STUDY,".pca_center.txt"),col.names=T,row.names=F,quote=F,sep="\t")

mylegend = c("Africa", "America", "East Asia", "Europe", "Middle-East", "South Asia")
# mycol = brewer.pal(6, "Dark2")
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E","#E7298A")
populations <- c("Africa", "American", "East Asia", "European", "Middle-East", "South Asia")
info$Color <- mycol[match(info$myREG, populations)]
#
x1=c(dataTGP$V3,mydata$V3); y1=c(dataTGP$V4,mydata$V4)
x2=c(dataTGP$V5,mydata$V5); y2=c(dataTGP$V6,mydata$V6)
x3=c(dataTGP$V6,mydata$V6); y3=c(dataTGP$V7,mydata$V7)
pdf(paste0(DATA_DIR,"/", STUDY,"/PCA_results/",STUDY,".pca.pdf"),width=8,height=4)
split.screen(c(1,2))
screen(1); plot(dataTGP$V3,dataTGP$V4 ,col=as.character(info$Color),pch=info$PCH,xlab="PC1",ylab="PC2",xlim=c(min(x1),max(x1)),ylim=c(min(y1),max(y1))); points(mydata$V3,mydata$V4,pch=16)
screen(2); plot(dataTGP$V5,dataTGP$V6 ,col=as.character(info$Color),pch=info$PCH,xlab="PC3",ylab="PC4",xlim=c(min(x2),max(x2)),ylim=c(min(y2),max(y2))); points(mydata$V5,mydata$V6,pch=16)
# screen(3); plot(dataTGP$V6,dataTGP$V7 ,col=as.character(info$Color),pch=info$PCH,xlab="PC4",ylab="PC5",xlim=c(min(x3),max(x3)),ylim=c(min(y3),max(y3))); points(mydata$V6,mydata$V7,pch=16)
# legend("topright",cex=0.5,legend=mylegend,col=mycol,pch=16)
dev.off()


#random forest
library(randomForest)
ancestry=as.character(info$myREG)
ancestry[which(ancestry=="Finnish")] = "European"
dataRF = cbind(ancestry,dataTGP[,3:tokeep])
model = randomForest(as.factor(ancestry) ~ ., data=dataRF)
# model = randomForest(as.factor(ancestry) ~ ., data=dataRF, ntree=500, mtry=4, nodesize=5)
pred  = predict(model, mydata[,3:tokeep], type="prob")
if (!dir.exists(paste0(DATA_DIR,"/", STUDY,"/rf_results"))) {
  dir.create(paste0(DATA_DIR,"/", STUDY,"/rf_results"))
}
library(dplyr)
pred_df <- as.data.frame(pred)
rf_ancestry <- pred_df %>%
  rowwise() %>%
  mutate(ancestry = names(.)[which.max(c_across(everything()))])
rownames(rf_ancestry) = rownames(pred)
write.table(rf_ancestry,file=paste0(DATA_DIR,"/", STUDY,"/rf_results/",STUDY,".pca_rf.txt"),col.names=T,row.names=T,quote=F,sep="\t")

