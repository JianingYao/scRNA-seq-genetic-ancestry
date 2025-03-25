#!/bin/bash

# scSNPs from a specific study
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for ERROR in $(seq 0.00 0.01 0.12); do
     for STUDY in "${list[@]}"; do
          Rscript plot_power.r $STUDY 'scSNPs' $ERROR
     done
done