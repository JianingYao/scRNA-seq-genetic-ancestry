#!/bin/bash

RESULT_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results"
TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
TGP_PLINK_HG38="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc"
INFO="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"
PLINK2="/project/gazal_569/soft/plink1.9/plink"

# Create exon_UTRS for QC
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

# QC the datasets: Finnish is considered as European
# grep -v oth $INFO | grep -v HGDP00621 | grep -v HGDP01302 | grep -v HGDP01261 | grep -v HGDP01270 | grep -v HGDP01271 | grep -v NA20314 | grep -v NA19625 | grep -v HG01108 | grep -v HG01241 | grep -v HG01242 | grep -v HG01243 | grep -v HG00277 | grep -v HG00344 | grep -v HG00350 | grep -v HG00380 | grep -v russian | grep -v uygur | grep -v hazara | grep -v clm | grep -vw pur | grep -vw mxl > info.txt
cut -f1-2 info.txt > $RESULT_DIR/idtokeep.txt
cut -f3-5 info.txt | sort | uniq | grep -v POP > $RESULT_DIR/listpop.txt

# Vary MAF and LD pruning threshold
for MAF in 05 02 01; do
  echo "===> MAF = 0.${MAF}"
  ############################
  # Build: hg19 / base build #
  ############################
  PREFIX="${RESULT_DIR}/TGP_PLINK_maf${MAF}"
  $PLINK2 --bfile "$TGP_PLINK" --maf 0.${MAF} --keep "$RESULT_DIR/idtokeep.txt" --make-bed --out "$PREFIX"
  # annotate base BIM + intersect
  PCA_FILE="$PREFIX"
  awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' "$PCA_FILE.bim" > "$PCA_FILE.bim.bed"
  bedtools intersect -a "$PCA_FILE.bim.bed" -b "$RESULT_DIR/exon_UTRS.bed" > "$PCA_FILE.exon_UTRS.bed"

  # LD pruning for multiple r^2 thresholds for allSNPs
  for R2 in 0.10 0.15 0.20; do
    TAG="r2_${R2/./p}"     # 0.10 -> 0p10, 0.15 -> 0p15, 0.20 -> 0p20
    OUTBASE="${PREFIX}.allSNPs_pruned.${TAG}"
    $PLINK2 --bfile "$PREFIX" --indep-pairwise 50 10 "$R2" --out "$OUTBASE"
    $PLINK2 --bfile "$PREFIX" --extract "${OUTBASE}.prune.in" --make-bed --out "$OUTBASE"
  done
  ############################
  # Build: hg38              #
  ############################
  PREFIX38="${RESULT_DIR}/TGP_PLINK_HG38_maf${MAF}"
  $PLINK2 --bfile "$TGP_PLINK_HG38" --maf 0.${MAF} --keep "$RESULT_DIR/idtokeep.txt" --make-bed --out "$PREFIX38"
  PCA_FILE="$PREFIX38"
  awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' "$PCA_FILE.bim" > "$PCA_FILE.bim.bed"
#   bedtools intersect -a "$PCA_FILE.bim.bed" -b "$RESULT_DIR/exon_UTRS.bed" > "$PCA_FILE.exon_UTRS.bed"

  echo "<=== done MAF = 0.${MAF}"
done


