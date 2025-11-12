#!/bin/bash

VCF_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/VCF"
list=("eQTLAutoimmune" "GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia" "SpatialMapSkin" "NSCL_lesions_tumor_classification" "nasalMucosaLifespan" "HumanCellLandscape" "GTEx_donor_tissue" "GTEx_donor" "HormoneRegulatedNetworksBreast" "humanHeartFailureCellularLandscape" "humanHeartFailureCellularLandscape_cell" "humanHeartFailureCellularLandscape_nuclei")

for STUDY in "${list[@]}"; do
  VCF="${VCF_DIR}/${STUDY}.SNP.Filtered.SV.vcf.gz"
  for MAF in 05 02 01; do                 
    for R2 in 0.10 0.15 0.20; do               
      R2_TAG="r2_${R2/./p}"                    
      DIR_NAME="Results_maf${MAF}_${R2_TAG}"
      if [ ! -d "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME" ]; then
        mkdir -p "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/$DIR_NAME"
      fi
      bash qc_pca_admixture.sh "$STUDY" "$VCF" "$MAF" "$R2" "$DIR_NAME"
    done
  done
done



