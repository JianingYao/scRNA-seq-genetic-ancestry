#!/bin/bash

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis"

# Power analysis using all common SNPs
allSNPs_DIR=$POWER_DIR/Results/allSNPs
if [ ! -d "$allSNPs_DIR" ]; then
    mkdir -p "$allSNPs_DIR"
fi
TGP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP"
cp $TGP_DIR/Results/*prune* $allSNPs_DIR
cd $allSNPs_DIR 
SNP="allSNPs"
while read LINE; do
    ARR=( $LINE ); POP=${ARR[0]};
    echo $POP
    sbatch --mem=30000 -t 0-03:00 --wrap="bash $POWER_DIR/pop_power.sh $SNP $POP"
done < $TGP_DIR/Results/listpop.txt

# Power analysis using sc SNPs from a HCA study
cd $POWER_DIR
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for STUDY in "${list[@]}"; do
    bash $POWER_DIR/scSNP_power.sh $STUDY $VCF
done

