#!/bin/bash

set -euo pipefail
shopt -s nullglob

POWER_DIR="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/05_latin/power_analysis"
INFO="/project2/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"

BASE_INFO="/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/info.txt"
IID_COL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if(tolower($i)=="iid"){print i; exit}}' "$BASE_INFO")

for POP in pur mxl clm; do
  awk -v IGNORECASE=1 -v TARGET="$POP" -F'\t' '
    NR==1 { for (i=1;i<=NF;i++) if (tolower($i)=="pop") pop=i; next }
    tolower($pop)==TARGET
  ' "$INFO" > "$POWER_DIR/$POP.txt"
  cat "$BASE_INFO" "$POWER_DIR/$POP.txt" > "$POWER_DIR/info_${POP}.tmp"
  {
    head -n1 "$POWER_DIR/info_${POP}.tmp"
    tail -n +2 "$POWER_DIR/info_${POP}.tmp" | sort -t $'\t' -k${IID_COL},${IID_COL}
  } > "$POWER_DIR/info_${POP}.txt"
  rm -f "$POWER_DIR/info_${POP}.tmp"
  cut -f1-2 $POWER_DIR/info_${POP}.txt > $POWER_DIR/idtokeep_$POP.txt
done

MAF=05
R2=0.10
R2_TAG="r2_${R2/./p}"
DIR_NAME="Results_maf${MAF}_${R2_TAG}"
if [ ! -d "$POWER_DIR/$DIR_NAME" ]; then
    mkdir -p "$POWER_DIR/$DIR_NAME"
fi

# Power analysis using all common SNPs
STUDY="allSNPs"
bash $POWER_DIR/allSNP_power.sh $STUDY $DIR_NAME $MAF $R2

# Power analysis using sc SNPs from a HCA study
cd "$POWER_DIR"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")
for STUDY in "${list[@]}"; do
    bash $POWER_DIR/scSNP_power.sh $STUDY $DIR_NAME $R2
done


