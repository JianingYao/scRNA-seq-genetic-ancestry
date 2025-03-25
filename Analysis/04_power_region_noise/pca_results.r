param <- commandArgs(trailingOnly=T)

ANC  = eval(paste(text=param[1])) 
SNPset = eval(paste(text=param[2]))
REGION = eval(paste(text=param[3]))
#ANC="European"; SNPset="allSNPs"; REGION="adygei"

info = read.table("../../myinfo.split.txt",h=T,sep="\t")
info = subset(info,info$myREG==ANC)

list_pop = unique(info$POP)

pcatokeep = c(7,5,5,5,7,7)
pcatokeep = rep(7,6)
names(pcatokeep) = c("Africa", "American", "EastAsia", "European", "MiddleEast", "SouthAsia")
tokeep = pcatokeep[ANC]

plotpca <-function(cluster,pos="topright"){
    mypop = unique(myinfo1$POP)
    #
    mycol1 = rep("grey",nrow(mydata1))
    for (i in 1:length(mypop)){
        mycol1[which(myinfo1$POP==mypop[i])] = i
    }
    mycol1 = as.numeric(mycol1)
    #
    mycol2 = rep("grey",nrow(mydata2))
    for (i in 1:length(mypop)){
        mycol2[which(myinfo2$POP==mypop[i])] = i
    }
    mycol2 = as.numeric(mycol2)
    #
    pdf(paste0("final_outs/",ANC,".", REGION, ".", cluster,".",SNPset,".pca.pdf"))
    plot  (mydata1$V3,mydata1$V4,col=mycol1,pch=16,xlab="PC1",ylab="PC2"); legend(pos,mypop,pch=16,col=1:length(mypop))
    points(mydata2$V3,mydata2$V4,col=mycol2,pch=1)
    plot  (mydata1$V4,mydata1$V5,col=mycol1,pch=16,xlab="PC2",ylab="PC3");
    points(mydata2$V4,mydata2$V5,col=mycol2,pch=1)
    plot  (mydata1$V5,mydata1$V6,col=mycol1,pch=16,xlab="PC3",ylab="PC4");
    points(mydata2$V5,mydata2$V6,col=mycol2,pch=1)
    plot  (mydata1$V6,mydata1$V7,col=mycol1,pch=16,xlab="PC4",ylab="PC5");
    points(mydata2$V6,mydata2$V7,col=mycol2,pch=1)
    dev.off()
}

#Cluster 1 analysis
num_clusters <- 10
for (cluster in 1:num_clusters){

    data=read.table(paste0("pca_results/",SNPset,".",ANC, ".", REGION, ".cluster",cluster,".eigenvec"),h=F)
    rownames(data) = data$V2
    
    cluster1 = which(info$NEWPOP != paste0(REGION, "_", cluster))
    cluster2 = which(info$NEWPOP == paste0(REGION, "_", cluster))
    mydata1  = data[cluster1,]
    myinfo1  = info[cluster1,]
    mydata2  = data[cluster2,]
    myinfo2  = info[cluster2,]

    #plot pca
    plotpca(cluster)

    #compute center of each reference population
    center = NULL
    for (pop in list_pop){
        center = rbind(center,
            apply(subset(mydata1[,3:tokeep],myinfo1$POP==pop),2,mean)
        )
    }
    rownames(center) = list_pop

    out=NULL
    thisiidpop=NULL
    for (i in 1:nrow(mydata2)){
        thisiid = NULL
        for (j in 1:nrow(center)){
           thisiid = c(thisiid,sum((mydata2[i,3:tokeep]-center[j,])**2))    
        }
        thisiidpop = c(thisiidpop,as.character(list_pop[which(thisiid == min(thisiid))]))
        thisiid = as.integer(thisiid == min(thisiid))
        out=rbind(out,thisiid)
    }
    colnames(out) = list_pop
    rownames(out) = mydata2[,2]
    table(as.character(myinfo2$POP),thisiidpop)
    write.table(table(as.character(myinfo2$POP),thisiidpop),file=paste0("pca_results/",ANC,".", REGION, ".", cluster,".",SNPset,".pca_center.txt"),col.names=T,row.names=T,quote=F,sep="\t")

    #
    library(randomForest)
    set.seed(100)
    ancestry = as.character(myinfo1$POP)
    dataRF = cbind(ancestry,mydata1[,3:tokeep])
    model <- randomForest(as.factor(ancestry) ~ ., data=dataRF)
    pred  <- predict(model, mydata2[,3:tokeep], type="prob")
    thisiidpop=NULL
    for (i in 1:nrow(pred)){
        thisiidpop = c(thisiidpop,colnames(pred)[head(which(pred[i,] == max(pred[i,])),1)])
    }
    mean.pop = NULL
    for (pop in colnames(pred)){
        mean.pop =rbind(mean.pop,
            apply(subset(pred,myinfo2$POP==pop),2,mean)
        )
    }
    out = table(as.character(myinfo2$POP),thisiidpop)
    out
    write.table(out,file=paste0("pcaRF_results/",ANC,".", REGION, ".", cluster,".",SNPset,".pca_randomforest.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}

