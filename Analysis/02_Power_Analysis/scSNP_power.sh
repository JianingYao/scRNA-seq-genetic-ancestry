#!/bin/bash

STUDY=$1
DIR_NAME=$2
R2=$3

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis"
STUDY_DIR=$POWER_DIR/$DIR_NAME/${STUDY}
if [ ! -d "$STUDY_DIR" ]; then
    mkdir -p "$STUDY_DIR"
fi
cd "$STUDY_DIR" 
SNP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME/${STUDY}"
TGP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP"
TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
PLINK2="/project/gazal_569/soft/plink1.9/plink"

#scSNPs analyses: Prune SNPs for PCA and Admixture
$PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $TGP_DIR/Results/idtokeep.txt --make-bed --out TGP_PLINK
$PLINK2 --bfile TGP_PLINK --indep-pairwise 50 10 $R2 --out TGP_PLINK.scSNPs_pruned
$PLINK2 --bfile TGP_PLINK --extract TGP_PLINK.scSNPs_pruned.prune.in --make-bed --out TGP_PLINK.scSNPs_pruned
#clean
rm $STUDY_DIR/TGP_PLINK.nosex $STUDY_DIR/TGP_PLINK.bed $STUDY_DIR/TGP_PLINK.bim $STUDY_DIR/TGP_PLINK.fam

# Run power analysis for each population
SNP="scSNPs"
while read LINE; do
    ARR=( $LINE ); POP=${ARR[0]};
    echo $POP
    sbatch --mem=30000 -t 1-00:00 -p main --wrap="bash $POWER_DIR/pop_power.sh $SNP $POP"
done < $TGP_DIR/Results/listpop.txt