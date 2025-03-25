param <- commandArgs(trailingOnly=T)

STUDY = eval(paste(text=param[1])) 

library(RColorBrewer)

DATA_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results/"
INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
data=read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".TGP_HGDP.pca.eigenvec"),h=F)
rownames(data) = data$V2
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
rownames(info) = info$IID
mydata  = data[as.character(read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".plink.fam"))[,2]),]
dataTGP = data[as.character(info$IID),]

mylegend = c("Africa", "America", "East Asia", "Europe", "Middle-East", "South Asia")
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E","#E7298A")
populations <- c("Africa", "American", "East Asia", "European", "Middle-East", "South Asia")
info$Color <- mycol[match(info$myREG, populations)]
#
x1=c(dataTGP$V3,mydata$V3); y1=c(dataTGP$V4,mydata$V4)
x2=c(dataTGP$V5,mydata$V5); y2=c(dataTGP$V6,mydata$V6)
x3=c(dataTGP$V6,mydata$V6); y3=c(dataTGP$V7,mydata$V7)
pdf(paste0(STUDY,".pca.pdf"),width=8,height=4)
split.screen(c(1,2))
screen(1); plot(dataTGP$V3,dataTGP$V4 ,col=as.character(info$Color),pch=info$PCH,xlab="PC1",ylab="PC2",xlim=c(min(x1),max(x1)),ylim=c(min(y1),max(y1))); points(mydata$V3,mydata$V4,pch=16)
screen(2); plot(dataTGP$V5,dataTGP$V6 ,col=as.character(info$Color),pch=info$PCH,xlab="PC3",ylab="PC4",xlim=c(min(x2),max(x2)),ylim=c(min(y2),max(y2))); points(mydata$V5,mydata$V6,pch=16)
dev.off()


