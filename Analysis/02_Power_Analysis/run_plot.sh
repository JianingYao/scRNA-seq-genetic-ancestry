#!/bin/bash

# all common SNPs
sbatch --mem=30000 -t 0-03:00 --wrap="Rscript plot_power.r 'allSNPs' 'allSNPs'"

# scSNPs from a specific study
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for STUDY in "${list[@]}"; do
     sbatch --mem=30000 -t 0-03:00 --wrap="Rscript plot_power.r $STUDY 'scSNPs'"
done