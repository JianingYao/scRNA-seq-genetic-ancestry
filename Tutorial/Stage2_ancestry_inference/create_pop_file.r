param <- commandArgs(trailingOnly=T)

STUDY = eval(paste(text=param[1])) 
DATA_DIR="../"
INFO_DIR="../../Analysis/00_TGP/"

bed=read.table(paste0(DATA_DIR,"/", STUDY,"/",STUDY,".HGDP_1kGP.fam"),h=F)[,1:2]
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
pop=rep(".",nrow(bed)); names(pop)=bed$V2

for (mypop in c("eas","amr","sas","afr","mid","eur")){
    pop[as.character(info$IID[which(info$REG==mypop)])] = mypop
}

write.table(pop,file=paste0(DATA_DIR,"/", STUDY,"/",STUDY,".HGDP_1kGP.pca.pop"),sep="\t",quote=F,col.names=F,row.names=F)

