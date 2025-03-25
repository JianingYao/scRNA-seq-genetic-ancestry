#!/bin/bash

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for ERROR in $(seq 0.00 0.01 0.12); do
    cd $POWER_DIR
    if [ ! -d "$POWER_DIR/Results_$ERROR" ]; then
    mkdir -p "$POWER_DIR/Results_$ERROR"
    fi
    for STUDY in "${list[@]}"; do
        STUDY_DIR=$POWER_DIR/Results_$ERROR/${STUDY}
        RESULT_DIR=$POWER_DIR/Results_$ERROR/${STUDY}/pca_results
        if [ ! -d "$STUDY_DIR" ]; then
            mkdir -p "$STUDY_DIR"
        fi
        if [ ! -d "$RESULT_DIR" ]; then
            mkdir -p "$RESULT_DIR"
        fi
        if [ ! -d "$STUDY_DIR/pcaRF_results" ]; then
            mkdir -p "$STUDY_DIR/pcaRF_results"
        fi
        if [ ! -d "$STUDY_DIR/admixture_results" ]; then
            mkdir -p "$STUDY_DIR/admixture_results"
        fi
        if [ ! -d "$STUDY_DIR/final_outs" ]; then
            mkdir -p "$STUDY_DIR/final_outs"
        fi
        cd $STUDY_DIR 
        for ANC in Africa EastAsia European MiddleEast SouthAsia American; do
            sbatch --mem=50000 -t 0-03:00 --partition=main --wrap="bash $POWER_DIR/power_region.sh $STUDY $ERROR $ANC"
        done
    done
done





