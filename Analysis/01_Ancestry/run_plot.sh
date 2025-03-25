#!/bin/bash

list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")

for STUDY in "${list[@]}"; do
    Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/pca_analysis.r $STUDY
    Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/plot_admixture.r $STUDY
done
