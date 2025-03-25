#!/bin/bash

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for ERROR in $(seq 0.00 0.01 0.12); do
    cd $POWER_DIR
    for STUDY in "${list[@]}"; do
        bash $POWER_DIR/scSNP_power.sh $STUDY $ERROR
    done
done

