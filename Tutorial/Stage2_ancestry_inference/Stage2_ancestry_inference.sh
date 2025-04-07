#!/bin/bash

STUDY=$1
VCF=$2

PLINK2="./plink1.9/plink"
# Note: Our recommended reference file is available at HGDP_1kGP.exon_UTRS.tar.xz. It contains 3,481 HGDP+1kGP individuals, and 302,006 common SNPs in coding exons and untranslated regions.
PCA_FILE="HGDP_1kGP.exon_UTRS"
tar -xJvf $PCA_FILE.tar.xz

# Step 1: Convert your vcf into plink format and QC the data
DIR="../${STUDY}"

if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"
fi

$PLINK2 --const-fid 0 --vcf $VCF --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step1
$PLINK2 --bfile $DIR/$STUDY.plink.step1 --mind 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step2
$PLINK2 --const-fid 0 --vcf $VCF --keep $DIR/$STUDY.plink.step2.fam --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink
rm $DIR/$STUDY.plink.step*

# Step 2: Find rs ID and create a list with rs id in both scRNA dataset and the reference dataset
perl annotate_rs_vcf.pl $DIR/$STUDY.plink $PCA_FILE
grep -v toremove $DIR/$STUDY.plink.bim | cut -f2 > $DIR/$STUDY.plink.list

# Step 3: Merge the 2 files and perform pruning
$PLINK2 --bfile $DIR/$STUDY.plink --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp1 --allow-extra-chr
$PLINK2 --bfile $PCA_FILE         --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp2
$PLINK2 --bfile $DIR/$STUDY.tmp1  --bmerge $DIR/$STUDY.tmp2 --make-bed --allow-no-sex --out $DIR/$STUDY.HGDP_1kGP
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP   --indep-pairwise 50 10 0.1 --out $DIR/$STUDY.HGDP_1kGP.pca
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP   --extract $DIR/$STUDY.HGDP_1kGP.pca.prune.in --make-bed --out $DIR/$STUDY.HGDP_1kGP.pca
rm $DIR/$STUDY.tmp*

# Step 4: Run PCA
awk '{print $1, $2, "MYID"}' $DIR/$STUDY.plink.fam > $DIR/$STUDY.HGDP_1kGP.cluster
awk '{print $1, $2, "HGDP_1kGP"}' $PCA_FILE.fam >> $DIR/$STUDY.HGDP_1kGP.cluster
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP.pca --pca --within $DIR/$STUDY.HGDP_1kGP.cluster --pca-cluster-names HGDP_1kGP --out $DIR/$STUDY.HGDP_1kGP.pca
Rscript pca_analysis.r $STUDY

# Step 5: Run ADMIXTURE
if [ ! -d "$DIR/ADMIXTURE_results" ]; then
    mkdir -p "$DIR/ADMIXTURE_results"
fi
Rscript create_pop_file.r $STUDY
cd $DIR/ADMIXTURE_results 
./bin/admixture_linux-1.3.0/admixture --supervised ../$DIR/$STUDY.HGDP_1kGP.pca.bed 6
Rscript ../../Stage2_ancestry_inference/admixture_analysis.r $STUDY
