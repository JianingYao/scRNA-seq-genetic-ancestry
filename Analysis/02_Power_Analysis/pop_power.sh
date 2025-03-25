#!/bin/bash

TGP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis"
INFO=$TGP_DIR/info.txt
PLINK2="/project/gazal_569/soft/plink1.9/plink"
module load r/4.3.2

SNP=$1
POP=$2

if [ ! -d "pca_results" ]; then
    mkdir -p "pca_results"
fi
if [ ! -d "pcaRF_results" ]; then
    mkdir -p "pcaRF_results"
fi
if [ ! -d "admixture_results" ]; then
    mkdir -p "admixture_results"
fi

#PCA
awk -v pop=$POP '(NR>1) {if($3==pop) {print $1, $2, "MYPOP"} else {print $1, $2, "TGP"}}' $INFO > TGP_PLINK.no$POP.cluster
$PLINK2 --bfile TGP_PLINK.${SNP}_pruned  --pca --within TGP_PLINK.no$POP.cluster --pca-cluster-names TGP --out TGP_PLINK.${SNP}_pca.no$POP
Rscript $POWER_DIR/pca_power.r $POP $SNP
#clean
rm TGP_PLINK.no$POP.cluster TGP_PLINK.*SNPs_pca.no$POP.eigen* TGP_PLINK.*SNPs_pca.no$POP.log TGP_PLINK.*SNPs_pca.no$POP.nosex

#Admixture
awk -v pop=$POP '(NR>1) {if($3==pop) {print "."} else {print $4}}' $INFO > TGP_PLINK.${SNP}_admixture.no$POP.pop
cp TGP_PLINK.${SNP}_pruned.bed  TGP_PLINK.${SNP}_admixture.no$POP.bed
cp TGP_PLINK.${SNP}_pruned.bim  TGP_PLINK.${SNP}_admixture.no$POP.bim
cp TGP_PLINK.${SNP}_pruned.fam  TGP_PLINK.${SNP}_admixture.no$POP.fam
/project/gazal_569/jianing/sc-RNA_ancestry/bin/admixture_linux-1.3.0/admixture --supervised TGP_PLINK.${SNP}_admixture.no$POP.bed 6
paste TGP_PLINK.${SNP}_admixture.no$POP.pop TGP_PLINK.${SNP}_admixture.no$POP.6.Q | awk '{if($1==".") {print $2, $3, $4, $5, $6, $7, $8} }' > admixture_results/${SNP}.no$POP.txt
#clean
rm TGP_PLINK.*SNPs_admixture.no$POP.*

