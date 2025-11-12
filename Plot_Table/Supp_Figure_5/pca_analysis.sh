#!/bin/bash

## without noise
PLINK2="/project/gazal_569/soft/plink1.9/plink"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region"

list=("allSNPs" "GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")
for STUDY in "${list[@]}"; do
    STUDY_DIR=$POWER_DIR/Results/${STUDY}
    RESULT_DIR=$POWER_DIR/Results/${STUDY}/pca_results
    for ANC in Africa EastAsia European MiddleEast SouthAsia American; do
        #PCA ALL
        $PLINK2 --bfile $STUDY_DIR/TGP_PLINK.$ANC.${STUDY}_pruned --pca --out $RESULT_DIR/${STUDY}.$ANC
    done
done


## with noise
PLINK2="/project/gazal_569/soft/plink1.9/plink"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise"

list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")
for STUDY in "${list[@]}"; do
    STUDY_DIR=$POWER_DIR/Results_0.08/${STUDY}
    RESULT_DIR=$POWER_DIR/Results_0.08/${STUDY}/pca_results
    for ANC in Africa EastAsia European MiddleEast SouthAsia American; do
        #PCA ALL
        $PLINK2 --bfile $STUDY_DIR/TGP_PLINK.$ANC.${STUDY}_pruned --pca --out $RESULT_DIR/${STUDY}.$ANC
    done
done
