param <- commandArgs(trailingOnly=T)
STUDY = eval(paste(text=param[1])) 
DIR_NAME = eval(paste(text=param[2])) 

DATA_DIR=paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/", DIR_NAME)
INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"

bed=read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".TGP.fam"),h=F)[,1:2]
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
pop=rep(".",nrow(bed)); names(pop)=bed$V2

for (mypop in c("eas","amr","sas","afr","mid","eur")){
    pop[as.character(info$IID[which(info$REG==mypop)])] = mypop
}

write.table(pop,file=paste0(DATA_DIR,"/", STUDY,"/",STUDY,".TGP_HGDP.pca.pop"),sep="\t",quote=F,col.names=F,row.names=F)

