#!/bin/bash

STUDY=$1
VCF=$2

PLINK2="./plink1.9/plink"
# Note: Follow the script at ../../Analysis/00_TGP/run_TGP.sh to first create and QC the reference HGDP+1kGP dataset, or select an alternative dataset of your choice.
PCA_FILE="../../Analysis/00_TGP/Results/TGP_PLINK"

# Step 1: Convert vcf into plink format and QC the data
DIR="./Ancestry/Results/${STUDY}"

if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"
fi

$PLINK2 --const-fid 0 --vcf $VCF --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step1
$PLINK2 --bfile $DIR/$STUDY.plink.step1 --mind 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step2
$PLINK2 --const-fid 0 --vcf $VCF --keep $DIR/$STUDY.plink.step2.fam --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink
rm $DIR/$STUDY.plink.step*

# Step 2: Find rs ID and create a list with rs id in both scRNA dataset and the reference dataset
perl annotate_rs_vcf.pl $DIR/$STUDY.plink $PCA_FILE.exon_UTRS
grep -v toremove $DIR/$STUDY.plink.bim | cut -f2 > $DIR/$STUDY.plink.list1
grep -v toremove $PCA_FILE.exon_UTRS.bed | cut -f4  > $DIR/$STUDY.plink.list2
join <(sort $DIR/$STUDY.plink.list1) <(sort $DIR/$STUDY.plink.list2) > $DIR/$STUDY.plink.list
rm $DIR/$STUDY.plink.list1 $DIR/$STUDY.plink.list2

# Step 3: Merge the 2 files and perform pruning
$PLINK2 --bfile $DIR/$STUDY.plink --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp1 --allow-extra-chr
$PLINK2 --bfile $PCA_FILE         --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp2
$PLINK2 --bfile $DIR/$STUDY.tmp1  --bmerge $DIR/$STUDY.tmp2 --make-bed --allow-no-sex --out $DIR/$STUDY.TGP
$PLINK2 --bfile $DIR/$STUDY.TGP   --indep-pairwise 50 10 0.1 --out $DIR/$STUDY.TGP_HGDP.pca
$PLINK2 --bfile $DIR/$STUDY.TGP   --extract $DIR/$STUDY.TGP_HGDP.pca.prune.in --make-bed --out $DIR/$STUDY.TGP_HGDP.pca
rm $DIR/$STUDY.tmp*

# Step 4: Run PCA
awk '{print $1, $2, "MYID"}' $DIR/$STUDY.plink.fam > $DIR/$STUDY.TGP_HGDP.cluster
awk '{print $1, $2, "TGP"}' $PCA_FILE.fam >> $DIR/$STUDY.TGP_HGDP.cluster
$PLINK2 --bfile $DIR/$STUDY.TGP_HGDP.pca --pca --within $DIR/$STUDY.TGP_HGDP.cluster --pca-cluster-names TGP --out $DIR/$STUDY.TGP_HGDP.pca
Rscript pca_analysis.r $STUDY

# Step 5: Run ADMIXTURE
if [ ! -d "$DIR/ADMIXTURE_results" ]; then
    mkdir -p "$DIR/ADMIXTURE_results"
fi
Rscript create_pop_file.r $STUDY
cd $DIR/ADMIXTURE_results 
./bin/admixture_linux-1.3.0/admixture --supervised $DIR/$STUDY.TGP_HGDP.pca.bed 6
Rscript admixture_analysis.r $STUDY
