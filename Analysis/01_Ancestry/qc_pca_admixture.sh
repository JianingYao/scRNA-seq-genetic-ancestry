#!/bin/bash

set -euo pipefail

STUDY=$1
VCF=$2
MAF=$3
R2=$4
DIR_NAME=$5

PLINK2="/project/gazal_569/soft/plink1.9/plink"
PCA_FILE="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results/TGP_PLINK_maf${MAF}"
DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME/${STUDY}"
if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"
fi

# step 1: convert vcf into plink format and QC the data
$PLINK2 --const-fid 0 --vcf $VCF --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step1
$PLINK2 --bfile $DIR/$STUDY.plink.step1 --mind 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step2
$PLINK2 --const-fid 0 --vcf $VCF --keep $DIR/$STUDY.plink.step2.fam --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink
rm $DIR/$STUDY.plink.step*

# step 2: find rs ID and create a list with rs id in both $FILE and in $PCA_FILE.exon_UTRS
perl /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/annotate_rs_vcf.pl $DIR/$STUDY.plink $PCA_FILE.exon_UTRS
grep -v toremove $DIR/$STUDY.plink.bim | cut -f2 > $DIR/$STUDY.plink.list1
grep -v toremove $PCA_FILE.exon_UTRS.bed | cut -f4  > $DIR/$STUDY.plink.list2
join <(sort $DIR/$STUDY.plink.list1) <(sort $DIR/$STUDY.plink.list2) > $DIR/$STUDY.plink.list
rm $DIR/$STUDY.plink.list1 $DIR/$STUDY.plink.list2

# step 3: merge the 2 files and do prunning
$PLINK2 --bfile $DIR/$STUDY.plink --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp1 --allow-extra-chr
$PLINK2 --bfile $PCA_FILE         --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp2
$PLINK2 --bfile $DIR/$STUDY.tmp1  --bmerge $DIR/$STUDY.tmp2 --make-bed --allow-no-sex --out $DIR/$STUDY.TGP
$PLINK2 --bfile $DIR/$STUDY.TGP   --indep-pairwise 50 10 $R2 --out $DIR/$STUDY.TGP_HGDP.pca
$PLINK2 --bfile $DIR/$STUDY.TGP   --extract $DIR/$STUDY.TGP_HGDP.pca.prune.in --make-bed --out $DIR/$STUDY.TGP_HGDP.pca
rm $DIR/$STUDY.tmp*

#run pca
awk '{print $1, $2, "MYID"}' $DIR/$STUDY.plink.fam > $DIR/$STUDY.TGP_HGDP.cluster
awk '{print $1, $2, "TGP"}' $PCA_FILE.fam >> $DIR/$STUDY.TGP_HGDP.cluster
$PLINK2 --bfile $DIR/$STUDY.TGP_HGDP.pca --pca --within $DIR/$STUDY.TGP_HGDP.cluster --pca-cluster-names TGP --out $DIR/$STUDY.TGP_HGDP.pca
Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/pca_analysis.r $STUDY $DIR_NAME

#run admixture
if [ ! -d "$DIR/ADMIXTURE_results" ]; then
    mkdir -p "$DIR/ADMIXTURE_results"
fi
Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/create_pop_file.r $STUDY $DIR_NAME
cd $DIR/ADMIXTURE_results 
/project/gazal_569/jianing/sc-RNA_ancestry/bin/admixture_linux-1.3.0/admixture --supervised $DIR/$STUDY.TGP_HGDP.pca.bed 6
Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/plot_admixture.r $STUDY $DIR_NAME
