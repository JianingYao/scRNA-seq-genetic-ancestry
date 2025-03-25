################ Step 1: Grouping populations ################ 
info = read.table("./Results/myinfo.QC0.txt",h=T,sep="\t")
info$POP = as.character(info$POP)

#ANC=="Africa"
    info$POP[which(info$POP=="lwk")]        = "kenya"
    info$POP[which(info$POP=="bantukenya")] = "kenya"
    info$POP[which(info$POP=="yri")]        = "nigeria_yoruba"
    info$POP[which(info$POP=="yoruba")]     = "nigeria_yoruba"
    info$POP[which(info$POP=="esn")]        = "nigeria_yoruba"
    info$POP[which(info$POP=="gwd")]        = "mandenka"
    info$POP[which(info$POP=="mandenka")]   = "mandenka" 
#ANC=="EastAsia"
    info$POP[which(info$POP=="cdx")]      = "dai_kinh"
    info$POP[which(info$POP=="khv")]      = "dai_kinh"
    info$POP[which(info$POP=="dai")]      = "dai_kinh" 
    info$POP[which(info$POP=="chb")]      = "china_han"
    info$POP[which(info$POP=="chs")]      = "china_han"
    info$POP[which(info$POP=="han")]      = "china_han"
    info$POP[which(info$POP=="japanese")] = "japan"
    info$POP[which(info$POP=="jpt")]      = "japan"
#ANC="European"
    info$POP[which(info$POP=="tsi")]      = "italy"
    info$POP[which(info$POP=="italian")]  = "italy"
    info$POP[which(info$POP=="tuscan")]   = "italy"
    info$POP[which(info$POP=="orcadian")] = "orcadian_gbr_ceu"
    info$POP[which(info$POP=="gbr")]      = "orcadian_gbr_ceu"
    info$POP[which(info$POP=="ceu")]      = "orcadian_gbr_ceu" 
#ANC="SouthAsia"
    info$POP[which(info$POP=="balochi")] = "balochi_brahui_makrani"
    info$POP[which(info$POP=="brahui")]  = "balochi_brahui_makrani"
    info$POP[which(info$POP=="makrani")] = "balochi_brahui_makrani"
    info$POP[which(info$POP=="pathan")]  = "burusho_pathan_sindhi"
    info$POP[which(info$POP=="sindhi")]  = "burusho_pathan_sindhi"
    info$POP[which(info$POP=="burusho")] = "burusho_pathan_sindhi"



plotpca <-function(ANC,pos="topright"){
    myinfo = subset(info,info$myREG==ANC)
    pca = read.table(paste0("./Results/allSNPs/pca_results/allSNPs0.",ANC,".eigenvec"),h=F)
    mypop = unique(myinfo$POP)
    #
    mycol = rep("grey",nrow(pca))
    for (i in 1:length(mypop)){
        mycol[which(myinfo$POP==mypop[i])] = i
    }
    mycol = as.numeric(mycol)
    pdf(paste0("./Results/allSNPs/pca_",ANC,".allSNPs.pdf"))
    plot(pca$V3,pca$V4,col=mycol,pch=16,xlab="PC1",ylab="PC2"); legend(pos,mypop,pch=16,col=1:length(mypop))
    plot(pca$V4,pca$V5,col=mycol,pch=16,xlab="PC2",ylab="PC3");
    plot(pca$V5,pca$V6,col=mycol,pch=16,xlab="PC3",ylab="PC4");
    plot(pca$V6,pca$V7,col=mycol,pch=16,xlab="PC4",ylab="PC5"); 
    dev.off()
}


plotpca("Africa", "bottomright")
plotpca("American")
plotpca("EastAsia")
plotpca("European")
plotpca("MiddleEast")
plotpca("SouthAsia","topleft")


################ Step 2: split each population in 2 ################
mytable=table(info$POP)
pop.to.keep = names(mytable)

info$NEWPOP = as.character(info$POP)
info$SPLIT  = 0

num_clusters <-10 
for (pop in pop.to.keep){
    pos = which(info$POP == pop)
    n = length(pos)
    set.seed(100)
    split = sample(1:n, n) 
    for (i in 1:num_clusters){
        idx_start = round((i-1) * n / num_clusters) + 1
        idx_end = round(i * n / num_clusters)
        info$NEWPOP[pos[split[idx_start:idx_end]]] = paste0(pop, "_", i)
        info$SPLIT[pos[split[idx_start:idx_end]]] = i
    }
}

write.table(info,file="./Results/myinfo.split.txt",col.names=T,row.names=F,sep="\t",quote=F)