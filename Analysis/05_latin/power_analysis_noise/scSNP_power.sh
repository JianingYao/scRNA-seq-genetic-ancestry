#!/bin/bash

STUDY=$1
ERROR=$2
DIR_NAME=$3
R2=$4

BASE_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/05_latin/power_analysis_noise"
POWER_DIR="${BASE_DIR}/${DIR_NAME}"
STUDY_DIR=$POWER_DIR/Results_$ERROR/${STUDY}
if [ ! -d "$POWER_DIR/Results_$ERROR" ]; then
    mkdir -p "$POWER_DIR/Results_$ERROR"
fi
if [ ! -d "$STUDY_DIR" ]; then
    mkdir -p "$STUDY_DIR"
fi
cd "$STUDY_DIR" 
SNP_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME/${STUDY}"
TGP_PLINK="/project2/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
PLINK2="/project2/gazal_569/soft/plink1.9/plink"

#scSNPs analyses: Prune SNPs for PCA and Admixture
for POP in pur mxl clm; do
    $PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $BASE_DIR/idtokeep_$POP.txt --make-bed --out TGP_PLINK_$POP
    $PLINK2 --bfile TGP_PLINK_$POP --indep-pairwise 50 10 $R2 --out TGP_PLINK_$POP.scSNPs_pruned
    $PLINK2 --bfile TGP_PLINK_$POP --extract TGP_PLINK_$POP.scSNPs_pruned.prune.in --make-bed --out TGP_PLINK_$POP.scSNPs_pruned
    #clean
    rm $STUDY_DIR/TGP_PLINK_$POP.nosex $STUDY_DIR/TGP_PLINK_$POP.bed $STUDY_DIR/TGP_PLINK_$POP.bim $STUDY_DIR/TGP_PLINK_$POP.fam

    # Run power analysis for each population
    SNP="scSNPs"
    echo $POP
    sbatch --mem=30000 -t 0-03:00 -p main --wrap="bash $BASE_DIR/pop_power.sh $SNP $POP $ERROR"
done
