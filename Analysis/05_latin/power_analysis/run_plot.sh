#!/bin/bash

set -euo pipefail
shopt -s nullglob

for MAF in 05 02 01; do
    for R2 in 0.10 0.15 0.20; do
        R2_TAG="r2_${R2/./p}"
        DIR_NAME="Results_maf${MAF}_${R2_TAG}"
        # all common SNPs
        sbatch --mem=30000 -t 1-00:00 --wrap="Rscript plot_power.r 'allSNPs' 'allSNPs' '$DIR_NAME'"

        # scSNPs from a specific study
        list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

        for STUDY in "${list[@]}"; do
          sbatch --mem=30000 -t 1-00:00 -p main --wrap="Rscript plot_power.r '$STUDY' 'scSNPs' '$DIR_NAME'"
        done
    done
done
