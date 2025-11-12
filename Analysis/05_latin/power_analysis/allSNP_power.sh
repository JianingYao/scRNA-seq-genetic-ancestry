#!/bin/bash

STUDY=$1 # allSNPs
DIR_NAME=$2
MAF=$3
R2=$4

POWER_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/05_latin/power_analysis"
TGP_PLINK="/project2/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
PLINK2="/project2/gazal_569/soft/plink1.9/plink"

STUDY_DIR=$POWER_DIR/$DIR_NAME/${STUDY}
if [ ! -d "$STUDY_DIR" ]; then
    mkdir -p "$STUDY_DIR"
fi
cd "$STUDY_DIR" 

for POP in pur mxl clm; do
    echo "===> MAF = 0.${MAF}"
    PREFIX="${STUDY_DIR}/TGP_PLINK_${POP}"
    $PLINK2 --bfile "$TGP_PLINK" --maf 0.${MAF} --keep "$POWER_DIR/idtokeep_$POP.txt" --make-bed --out "$PREFIX"
    # annotate base BIM + intersect
    PCA_FILE="$PREFIX"
    awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' "$PCA_FILE.bim" > "$PCA_FILE.bim.bed"
    bedtools intersect -a "$PCA_FILE.bim.bed" -b "/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results/exon_UTRS.bed" > "$PCA_FILE.exon_UTRS.bed"

    # LD pruning for multiple r^2 thresholds for allSNPs
    OUTBASE="${PREFIX}.allSNPs_pruned"
    $PLINK2 --bfile "$PREFIX" --indep-pairwise 50 10 "$R2" --out "$OUTBASE"
    $PLINK2 --bfile "$PREFIX" --extract "${OUTBASE}.prune.in" --make-bed --out "$OUTBASE"
    rm $STUDY_DIR/TGP_PLINK_${POP}.nosex $STUDY_DIR/TGP_PLINK_${POP}.bed $STUDY_DIR/TGP_PLINK_${POP}.bim $STUDY_DIR/TGP_PLINK_${POP}.fam

    # Run power analysis for each population
    SNP="allSNPs"
    echo $POP
    sbatch --mem=30000 -t 0-03:00 -p main --wrap="bash $POWER_DIR/pop_power.sh $SNP $POP"
done
