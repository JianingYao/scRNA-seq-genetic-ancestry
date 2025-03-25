#!/bin/bash

TGP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise"
INFO=$POWER_DIR/sorted_info.txt
# INFO=$POWER_DIR/info.txt
PLINK2="/project/gazal_569/soft/plink1.9/plink"
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.3.2

SNP=$1
POP=$2
ERROR=$3

if [ ! -d "pca_results" ]; then
    mkdir -p "pca_results"
fi
if [ ! -d "pcaRF_results" ]; then
    mkdir -p "pcaRF_results"
fi
if [ ! -d "admixture_results" ]; then
    mkdir -p "admixture_results"
fi

#create new dataset with noise
awk -v pop=$POP '(NR>1) {if($3==pop) {print $1, $2}}' $INFO > $POP.$SNP.list
$PLINK2 --bfile TGP_PLINK.${SNP}_pruned --remove $POP.$SNP.list --make-bed --out $POP.$SNP.tmp1
$PLINK2 --bfile TGP_PLINK.${SNP}_pruned --keep   $POP.$SNP.list --recode transpose --out $POP.$SNP.tmp2
Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise/add_noise_seed.r $POP.$SNP $ERROR
$PLINK2 --tfile $POP.$SNP.tmp2 --make-bed --out $POP.$SNP.tmp2
# $PLINK2 --bfile $POP.$SNP.tmp1 --bmerge $POP.$SNP.tmp2 --make-bed --out $POP.$SNP.tmp
$PLINK2 --bfile $POP.$SNP.tmp1 --bmerge $POP.$SNP.tmp2 --make-bed --indiv-sort f $POWER_DIR/sorted_idtokeep.txt --out $POP.$SNP.tmp
rm $POP.$SNP.tmp?.*

#PCA
awk -v pop=$POP '(NR>1) {if($3==pop) {print $1, $2, "MYPOP"} else {print $1, $2, "TGP"}}' $INFO > TGP_PLINK.no$POP.cluster
$PLINK2 --bfile $POP.$SNP.tmp --pca --within TGP_PLINK.no$POP.cluster --pca-cluster-names TGP --out TGP_PLINK.${SNP}_pca.no$POP
Rscript $POWER_DIR/pca_power.r $POP $SNP
#clean
rm TGP_PLINK.no$POP.cluster TGP_PLINK.*SNPs_pca.no$POP.eigen* TGP_PLINK.*SNPs_pca.no$POP.log TGP_PLINK.*SNPs_pca.no$POP.nosex

#Admixture
awk -v pop=$POP '(NR>1) {if($3==pop) {print "."} else {print $4}}' $INFO > TGP_PLINK.${SNP}_admixture.no$POP.pop
cp $POP.$SNP.tmp.bed  TGP_PLINK.${SNP}_admixture.no$POP.bed
cp $POP.$SNP.tmp.bim  TGP_PLINK.${SNP}_admixture.no$POP.bim
cp $POP.$SNP.tmp.fam  TGP_PLINK.${SNP}_admixture.no$POP.fam
/project/gazal_569/jianing/sc-RNA_ancestry/bin/admixture_linux-1.3.0/admixture --supervised TGP_PLINK.${SNP}_admixture.no$POP.bed 6
paste TGP_PLINK.${SNP}_admixture.no$POP.pop TGP_PLINK.${SNP}_admixture.no$POP.6.Q | awk '{if($1==".") {print $2, $3, $4, $5, $6, $7, $8} }' > admixture_results/${SNP}.no$POP.txt
#clean
rm TGP_PLINK.*SNPs_admixture.no$POP.* $POP.$SNP.*

