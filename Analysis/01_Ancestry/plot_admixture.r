param <- commandArgs(trailingOnly=T)

STUDY = eval(paste(text=param[1])) 

library(RColorBrewer)

INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")

DATA_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results/"
ADM_DIR=paste0(DATA_DIR, STUDY, "/ADMIXTURE_results/")
admixture = read.table(paste0(ADM_DIR, STUDY,".TGP_HGDP.pca.6.Q"),h=F)
admixture.pop = read.table(paste0(DATA_DIR, STUDY,"/",STUDY,".TGP_HGDP.pca.pop"))
admixture.ind = read.table(paste0(DATA_DIR, STUDY,"/",STUDY,".TGP_HGDP.pca.fam"))[,2]

colnames(admixture)=c("eur","eas","amr","sas","afr","mid")
rownames(admixture)=admixture.ind

myadmixture  = admixture[as.character(read.table(paste0(DATA_DIR, STUDY,"/",STUDY,".plink.fam"))[,2]),] # subset only mydata
toplot=myadmixture[,order(apply(myadmixture,2,mean),decreasing=T)]
mycol=c("#E6AB02","#7570B3","#D95F02","#E7298A","#1B9E77","#66A61E")
mycol=mycol[order(apply(myadmixture,2,mean),decreasing=T)]

for (i in 6:1){
    toplot=toplot[order(toplot[,i],decreasing=T),]
}
pdf(paste0(ADM_DIR, "/", STUDY,".admixture.pdf"),width=10,height=4)
barplot(t(as.matrix(toplot)), col=mycol,xlab="", ylab="Admixture proportion", border=NA,xaxt="n")
dev.off()







