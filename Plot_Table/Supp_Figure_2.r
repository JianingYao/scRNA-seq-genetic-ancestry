mylist=c("allSNPs","GompertsAirwatCfCells","HnsccImmuneLandscape","humanPreimplantationEmbryos","StemCellsInChronicMyeloidLeukemia","SpatialMapSkin","NSCL_lesions_tumor_classification","nasalMucosaLifespan","HumanCellLandscape","HormoneRegulatedNetworksBreast","humanHeartFailureCellularLandscape_cell","humanHeartFailureCellularLandscape_nuclei")
mynames=c("All SNPs", "Airway epithelium", "CD45+", "Embryo", "Bone marrow", "Skin", "Lung", "Nose", " Cell Landscape", "Breast", "Heart (cell)", "Heart (nuclei)")
 
data = NULL
for (i in mylist){
    data = cbind(data,read.table(paste0("/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/06_inform_snp/inf_analysis/", i,".inf.txt"))[,1])
}
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E")

pdf("Supp_Figure_2.pdf",height=6,width=10)
par(mar=c(12,4,4,1))
barplot(data,beside=T,log="y",ylim=c(100,500000),border=0,col = mycol,ylab="# ancestry informative SNPs",yaxt="n",names.arg=mynames,las=2)
axis(2, at=c(100, 1000, 10000, 9e4), labels=c("100","1,000","10,000","100,000"))
legend("topright",c("PC1","PC2","PC3","PC4","PC5"),fill=mycol,bty="n",border=0)
dev.off()
