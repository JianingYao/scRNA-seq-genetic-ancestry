param <- commandArgs(trailingOnly=T)

STUDY  = eval(paste(text=param[1])) 
ERROR = as.character(eval(paste(text=param[2]))) 

num_clusters = 10
info = read.table("./myinfo.split.txt",h=T)

#pca_randomforest
pca_randomforest.anc.accuracy.STUDY = NULL
pca_randomforest.pop.accuracy.STUDY = NULL
list_anc = c("Africa", "EastAsia", "European", "MiddleEast", "SouthAsia", "American")
for (anc in list_anc){
    anc_info = subset(info,info$myREG==anc)
    list_pop = sort(unique(anc_info$POP))
    data.full = matrix(0,length(list_pop),length(list_pop))
    rownames(data.full) = list_pop
    colnames(data.full) = list_pop
    for (REGION in list_pop){
        for (i in 1:num_clusters){
            data = read.table(paste0("./Results_", ERROR, "/",STUDY,"/pcaRF_results/",anc,".",REGION,".",i,".",STUDY,".pca_randomforest.txt"),h=T)
            data.full[rownames(data),colnames(data)] = as.matrix(data.full[rownames(data),colnames(data)]+data[rownames(data),colnames(data)])
        }
    }
    pca_randomforest.anc.accuracy.STUDY = c(pca_randomforest.anc.accuracy.STUDY,sum(diag(data.full))/sum(data.full))
    pca_randomforest.pop.accuracy.STUDY = c(pca_randomforest.pop.accuracy.STUDY,diag(data.full)/apply(data.full,1,sum))
}
names(pca_randomforest.anc.accuracy.STUDY) = list_anc
pca_randomforest.anc.accuracy.STUDY
pca_randomforest.pop.accuracy.STUDY
write.table(pca_randomforest.anc.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_RF_anc_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")
write.table(pca_randomforest.pop.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_RF_pop_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")

#admixture
myinfo = read.table("./myinfo.split.txt",h=T)
admixture.anc.accuracy.STUDY = NULL
admixture.pop.accuracy.STUDY = NULL
list_anc = c("Africa", "EastAsia", "European", "MiddleEast", "SouthAsia", "American")
for (anc in list_anc){
    anc_info = subset(info,info$myREG==anc)
    list_pop = sort(unique(anc_info$POP))
    data = read.table(paste0("./Results_", ERROR, "/",STUDY,"/admixture_results/",STUDY,".",anc,".",list_pop[1],".1.txt"),h=F)
    data.full = matrix(0,ncol(data),ncol(data))
    rownames(data.full) = unique(subset(myinfo$POP,myinfo$myREG==anc & myinfo$SPLIT!=1))
    colnames(data.full) = rownames(data.full)
    for (REGION in list_pop){
        for (i in 1:num_clusters){
            data = read.table(paste0("./Results_", ERROR, "/",STUDY,"/admixture_results/",STUDY,".",anc,".",REGION,".",i,".txt"),h=F)
            colnames(data) = unique(subset(myinfo$POP,myinfo$myREG==anc & !(myinfo$POP == REGION & myinfo$SPLIT == i)))
            data = data[,colnames(data.full)]
            true.pop = subset(myinfo$POP,myinfo$myREG==anc & myinfo$POP == REGION & myinfo$SPLIT==i)
            for (mypop in colnames(data.full)){
                for (j in 1:ncol(data.full)){
                    data.full[mypop,j] = data.full[mypop,j] + sum(data[which(true.pop == mypop),j] == apply(data[which(true.pop == mypop),],1,max))
                }
            }
        }
    }
    admixture.anc.accuracy.STUDY = c(admixture.anc.accuracy.STUDY,sum(diag(data.full))/sum(data.full))
    admixture.pop.accuracy.STUDY = c(admixture.pop.accuracy.STUDY,diag(data.full)/apply(data.full,1,sum))
}
names(admixture.anc.accuracy.STUDY) = list_anc
admixture.anc.accuracy.STUDY
admixture.pop.accuracy.STUDY
write.table(admixture.anc.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_admixture_anc_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")
write.table(admixture.pop.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_admixture_pop_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")


# pca_center
anc.accuracy.STUDY = NULL
pop.accuracy.STUDY = NULL
list_anc = c("Africa", "EastAsia", "European", "MiddleEast", "SouthAsia", "American")
for (anc in list_anc){
    anc_info = subset(info,info$myREG==anc)
    list_pop = sort(unique(anc_info$POP))
    data.full = matrix(0,length(list_pop),length(list_pop))
    rownames(data.full) = list_pop
    colnames(data.full) = list_pop
    for (REGION in list_pop){
        for (i in 1:num_clusters){
            data = read.table(paste0("./Results_", ERROR, "/",STUDY,"/pca_results/",anc,".",REGION,".",i,".",STUDY,".pca_center.txt"),h=T)
            data.full[rownames(data),colnames(data)] = as.matrix(data.full[rownames(data),colnames(data)]+data[rownames(data),colnames(data)])
        }
    }
    anc.accuracy.STUDY = c(anc.accuracy.STUDY,sum(diag(data.full))/sum(data.full))
    pop.accuracy.STUDY = c(pop.accuracy.STUDY,diag(data.full)/apply(data.full,1,sum))
}
names(anc.accuracy.STUDY) = list_anc
anc.accuracy.STUDY
pop.accuracy.STUDY
write.table(anc.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_pca_anc_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")
write.table(pop.accuracy.STUDY,file=paste0("./Results_", ERROR, "/", STUDY, "/final_outs/", STUDY, "_pca_pop_accuracy.txt"),col.names=T,row.names=T,quote=F,sep="\t")


