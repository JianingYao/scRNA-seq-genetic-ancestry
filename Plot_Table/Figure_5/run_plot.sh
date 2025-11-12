#!/bin/bash

list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

for STUDY in "${list[@]}"; do
    for MAF in 05; do
        for R2 in 0.10; do
            R2_TAG="r2_${R2/./p}"
            DIR_NAME="Results_maf${MAF}_${R2_TAG}"
            Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/pca_plot.r $STUDY $DIR_NAME
            Rscript /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/plot_admixture.r $STUDY
        done
    done
done
