#!/bin/bash

# scSNPs from a specific study
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

for MAF in 05 02 01; do
  echo "[$(date '+%F %T')] >> Start MAF=${MAF}"

  for R2 in 0.10 0.15 0.20; do
    R2_TAG="r2_${R2/./p}"
    DIR_NAME="Results_maf${MAF}_${R2_TAG}"
    echo "[$(date '+%F %T')]  > Start R2=${R2} (${R2_TAG}), DIR=${DIR_NAME}"

    for ERROR in $(seq -w 0.00 0.01 0.00); do
      echo "[$(date '+%F %T')]   - Start ERROR=${ERROR}"

      for STUDY in "${list[@]}"; do
        Rscript plot_power.r "$STUDY" "scSNPs" "$ERROR" "$DIR_NAME"
        echo "[$(date '+%F %T')]     done STUDY=${STUDY} (ERROR=${ERROR}, DIR=${DIR_NAME})"
      done

      echo "[$(date '+%F %T')]   Completed ERROR=${ERROR} for R2=${R2}, MAF=${MAF}"
    done

    echo "[$(date '+%F %T')]  Completed R2=${R2} (${R2_TAG}) for MAF=${MAF}"
  done

  echo "[$(date '+%F %T')] Completed MAF=${MAF}"
done

echo "[$(date '+%F %T')] All done."

