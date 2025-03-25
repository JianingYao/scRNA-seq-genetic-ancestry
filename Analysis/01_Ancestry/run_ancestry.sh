#!/bin/bash

VCF_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/VCF"
list=("GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")

for STUDY in "${list[@]}"; do
    VCF="${VCF_DIR}/${STUDY}.SNP.Filtered.SV.vcf.gz"
    bash qc_pca_admixture.sh $STUDY $VCF
done


