#!/bin/bash

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise"

list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for ERROR in $(seq 0.00 0.01 0.12); do
    for STUDY in "${list[@]}"; do
        Rscript $POWER_DIR/check_output.r $STUDY $ERROR
    done
    Rscript $POWER_DIR/plot_error.r $ERROR
done

