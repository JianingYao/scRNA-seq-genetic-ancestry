#!/bin/bash

STUDY=$1
ERROR=$2
DIR_NAME=$3
R2=$4

BASE_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise"
POWER_DIR="${BASE_DIR}/${DIR_NAME}"
STUDY_DIR=$POWER_DIR/Results_$ERROR/${STUDY}
if [ ! -d "$POWER_DIR/Results_$ERROR" ]; then
    mkdir -p "$POWER_DIR/Results_$ERROR"
fi
if [ ! -d "$STUDY_DIR" ]; then
    mkdir -p "$STUDY_DIR"
fi
cd "$STUDY_DIR" 
SNP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME/${STUDY}"
TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
PLINK2="/project/gazal_569/soft/plink1.9/plink"

#scSNPs analyses: Prune SNPs for PCA and Admixture
$PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $BASE_DIR/sorted_idtokeep.txt --make-bed --indiv-sort f $BASE_DIR/sorted_idtokeep.txt --out TGP_PLINK
$PLINK2 --bfile TGP_PLINK --indep-pairwise 50 10 $R2 --out TGP_PLINK.scSNPs_pruned
$PLINK2 --bfile TGP_PLINK --extract TGP_PLINK.scSNPs_pruned.prune.in --make-bed --out TGP_PLINK.scSNPs_pruned
#clean
rm $STUDY_DIR/TGP_PLINK.nosex $STUDY_DIR/TGP_PLINK.bed $STUDY_DIR/TGP_PLINK.bim $STUDY_DIR/TGP_PLINK.fam

# Run power analysis for each population
SNP="scSNPs"
while read LINE; do
    ARR=( $LINE ); POP=${ARR[0]};
    echo $POP
    sbatch --mem=30000 -t 1-00:00 -p main --wrap="bash $BASE_DIR/pop_power.sh $SNP $POP $ERROR"
done < $BASE_DIR/sorted_listpop.txt
