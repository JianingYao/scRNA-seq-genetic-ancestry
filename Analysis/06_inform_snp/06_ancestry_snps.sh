#!/bin/bash
module load plink2
TGP_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP"
TGP_PLINK="/project2/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
PLINK2="/project2/gazal_569/soft/plink1.9/plink"
INFO=$TGP_DIR/info.txt
cut -f1,2,4 $INFO > ancestryfile.txt #keep the column with the ancestry label
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

# allSNPs
$PLINK2 --bfile $TGP_PLINK --maf 0.05 --keep $TGP_DIR/Results/idtokeep.txt --make-bed --out TGP_PLINK
$PLINK2 --bfile TGP_PLINK --indep-pairwise 50 10 0.1 --out TGP_PLINK.allSNPs_pruned
$PLINK2 --bfile TGP_PLINK --extract TGP_PLINK.allSNPs_pruned.prune.in --make-bed --out TGP_PLINK.allSNPs_pruned
$PLINK2 --bfile TGP_PLINK.allSNPs_pruned --pca --out allSNPs
awk '{print $1, $2, $3, $4, $5, $6, $7}' allSNPs.eigenvec > pheno.txt
#assoc
plink2 --bfile TGP_PLINK.allSNPs_pruned --glm allow-no-covars --pheno pheno.txt --out allSNPs
for PC in {1..5}; do
    awk '{if ($12<5*10**-8){print $0}}' allSNPs.PHENO$PC.glm.linear | wc -l >> allSNPs.inf.txt
done

# scSNPs
for STUDY in "${list[@]}"; do
    SNP_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results_maf05_r2_0p10/${STUDY}"
    $PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $TGP_DIR/Results/idtokeep.txt --make-bed --out TGP_PLINK
    $PLINK2 --bfile TGP_PLINK --indep-pairwise 50 10 0.1 --out TGP_PLINK.scSNPs_pruned
    $PLINK2 --bfile TGP_PLINK --extract TGP_PLINK.scSNPs_pruned.prune.in --make-bed --out TGP_PLINK.scSNPs_pruned
    plink2 --bfile TGP_PLINK.scSNPs_pruned --glm allow-no-covars --pheno pheno.txt --out $STUDY
    for PC in {1..5}; do
        awk '{if ($12<5*10**-8){print $0}}' $STUDY.PHENO$PC.glm.linear | wc -l >> $STUDY.inf.txt
    done
done

R

mylist=c("allSNPs","GompertsAirwatCfCells","HnsccImmuneLandscape","humanPreimplantationEmbryos","StemCellsInChronicMyeloidLeukemia","SpatialMapSkin","NSCL_lesions_tumor_classification","nasalMucosaLifespan","HumanCellLandscape","HormoneRegulatedNetworksBreast","humanHeartFailureCellularLandscape_cell","humanHeartFailureCellularLandscape_nuclei")
mynames=c("All SNPs", "Airway epithelium", "CD45+", "Embryo", "Bone marrow", "Skin", "Lung", "Nose", " Cell Landscape", "Breast", "Heart (cell)", "Heart (nuclei)")

data = NULL
for (i in mylist){
    data = cbind(data,read.table(paste0(i,".inf.txt"))[,1])
}
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E")

pdf("inf.pdf",height=6,width=10)
par(mar=c(12,4,4,1))
barplot(data,beside=T,log="y",ylim=c(100,500000),border=0,col = mycol,ylab="# ancestry informative SNPs",yaxt="n",names.arg=mynames,las=2)
axis(2,at=c(100,1000,10000,100000),c("100","1,000","10,000","100,000"))
legend("topright",c("PC1","PC2","PC3","PC4","PC5"),fill=mycol,bty="n",border=0)
dev.off()

