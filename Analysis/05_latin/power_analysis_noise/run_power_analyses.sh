#!/bin/bash

set -euo pipefail
shopt -s nullglob

POWER_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/05_latin/power_analysis_noise"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

MAF=05
R2=0.10
R2_TAG="r2_${R2/./p}"
DIR_NAME="Results_maf${MAF}_${R2_TAG}"
if [ ! -d "$POWER_DIR/$DIR_NAME" ]; then
    mkdir -p "$POWER_DIR/$DIR_NAME"
fi
ERROR=0.08
cd "$POWER_DIR/$DIR_NAME"
for STUDY in "${list[@]}"; do
    bash $POWER_DIR/scSNP_power.sh $STUDY $ERROR $DIR_NAME $R2
done



