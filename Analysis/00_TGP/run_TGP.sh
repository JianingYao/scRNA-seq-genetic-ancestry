#!/bin/bash

RESULT_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results"
TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
INFO="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"
PLINK2="/project/gazal_569/soft/plink1.9/plink"

#QC the datasets: Finnish is considered as European
# grep -v oth $INFO | grep -v HGDP00621 | grep -v HGDP01302 | grep -v HGDP01261 | grep -v HGDP01270 | grep -v HGDP01271 | grep -v NA20314 | grep -v NA19625 | grep -v HG01108 | grep -v HG01241 | grep -v HG01242 | grep -v HG01243 | grep -v HG00277 | grep -v HG00344 | grep -v HG00350 | grep -v HG00380 | grep -v russian | grep -v uygur | grep -v hazara | grep -v clm | grep -vw pur | grep -vw mxl > info.txt
cut -f1-2 info.txt > $RESULT_DIR/idtokeep.txt
cut -f3-5 info.txt | sort | uniq | grep -v POP > $RESULT_DIR/listpop.txt

#allSNPs analyses: Prune SNPs for PCA and Admixture
$PLINK2 --bfile $TGP_PLINK --maf 0.05 --keep $RESULT_DIR/idtokeep.txt --make-bed --out $RESULT_DIR/TGP_PLINK
$PLINK2 --bfile $RESULT_DIR/TGP_PLINK --indep-pairwise 50 10 0.1 --out $RESULT_DIR/TGP_PLINK.allSNPs_pruned
$PLINK2 --bfile $RESULT_DIR/TGP_PLINK --extract $RESULT_DIR/TGP_PLINK.allSNPs_pruned.prune.in --make-bed --out $RESULT_DIR/TGP_PLINK.allSNPs_pruned
#clean
# rm *log *prune* *nosex


#create exon_UTRS for QC
QCDIR="/project/gazal_569/DATA/ldsc/reference_files/bed"
GENOME="/project/gazal_569/DATA/gene_info/hg19/genome_big"
cat $QCDIR/Coding_UCSC.bed > $RESULT_DIR/exon_UTRS.tmp1.bed
cat $QCDIR/UTR_3_UCSC.bed >> $RESULT_DIR/exon_UTRS.tmp1.bed
cat $QCDIR/UTR_5_UCSC.bed >> $RESULT_DIR/exon_UTRS.tmp1.bed
sort -k 1,1 -k2,2n $RESULT_DIR/exon_UTRS.tmp1.bed > $RESULT_DIR/exon_UTRS.tmp2.bed
module load bedtools2
bedtools slop -i $RESULT_DIR/exon_UTRS.tmp2.bed -g genome_big -b 150 > $RESULT_DIR/exon_UTRS.tmp3.bed
bedtools merge -i $RESULT_DIR/exon_UTRS.tmp3.bed > $RESULT_DIR/exon_UTRS.bed
rm $RESULT_DIR/exon_UTRS.tmp*

#annotate PCA_FILE
PCA_FILE="${RESULT_DIR}/TGP_PLINK"
awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' $PCA_FILE.bim > $PCA_FILE.bim.bed
bedtools intersect -a $PCA_FILE.bim.bed -b $RESULT_DIR/exon_UTRS.bed > $PCA_FILE.exon_UTRS.bed


