#!/bin/bash

set -euo pipefail
shopt -s nullglob

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis"

for MAF in 05 02 01; do
    for R2 in 0.10 0.15 0.20; do
        R2_TAG="r2_${R2/./p}"
        DIR_NAME="Results_maf${MAF}_${R2_TAG}"
        if [ ! -d "$POWER_DIR/$DIR_NAME" ]; then
            mkdir -p "$POWER_DIR/$DIR_NAME"
        fi

        # Power analysis using all common SNPs
        allSNPs_DIR=$POWER_DIR/$DIR_NAME/allSNPs
        if [ ! -d "$allSNPs_DIR" ]; then
            mkdir -p "$allSNPs_DIR"
        fi
        TGP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP"
        src_prefix="$TGP_DIR/Results/TGP_PLINK_maf${MAF}.allSNPs_pruned.${R2_TAG}"
        for src in "${src_prefix}"*; do
            base="${src##*/}"
            # keep the suffix after the original prefix
            suffix="${base#TGP_PLINK_maf${MAF}.allSNPs_pruned.${R2_TAG}}"
            cp -v "$src" "$allSNPs_DIR/TGP_PLINK.allSNPs_pruned${suffix}"
        done
        cd "$allSNPs_DIR" 
        SNP="allSNPs"
        while read LINE; do
            ARR=( $LINE ); POP=${ARR[0]};
            echo $POP
            sbatch --mem=30000 -t 1-00:00 --wrap="bash $POWER_DIR/pop_power.sh $SNP $POP"
        done < $TGP_DIR/Results/listpop.txt

        # Power analysis using sc SNPs from a HCA study
        cd "$POWER_DIR"
        list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")
        for STUDY in "${list[@]}"; do
            bash $POWER_DIR/scSNP_power.sh $STUDY $DIR_NAME $R2
        done
    done
done

