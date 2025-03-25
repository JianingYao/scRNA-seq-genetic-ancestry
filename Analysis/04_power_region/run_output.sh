#!/bin/bash

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region"

list=("allSNPs" "GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for STUDY in "${list[@]}"; do
    Rscript $POWER_DIR/check_output.r $STUDY
done

Rscript $POWER_DIR/plot_error.r