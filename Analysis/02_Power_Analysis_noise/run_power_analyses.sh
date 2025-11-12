#!/bin/bash

set -euo pipefail
shopt -s nullglob

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

for MAF in 05 02 01; do
    for R2 in 0.10 0.15 0.20; do
        R2_TAG="r2_${R2/./p}"
        DIR_NAME="Results_maf${MAF}_${R2_TAG}"
        if [ ! -d "$POWER_DIR/$DIR_NAME" ]; then
            mkdir -p "$POWER_DIR/$DIR_NAME"
        fi
        for ERROR in $(seq 0.00 0.01 0.13); do
            cd "$POWER_DIR/$DIR_NAME"
            for STUDY in "${list[@]}"; do
                bash $POWER_DIR/scSNP_power.sh $STUDY $ERROR $DIR_NAME $R2
            done
        done
    done
done

