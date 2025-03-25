#!/bin/bash

list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")

for STUDY in "${list[@]}"; do
    Rscript ../Figure_3/pca_analysis.r $STUDY
    Rscript ../Figure_3/plot_admixture.r $STUDY
done
