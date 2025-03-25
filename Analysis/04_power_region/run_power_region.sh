#!/bin/bash

TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
INFO="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"
PLINK2="/project/gazal_569/soft/plink1.9/plink"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region"
if [ ! -d "$POWER_DIR/Results" ]; then
    mkdir -p "$POWER_DIR/Results"
fi
# module load r/4.3.2

################ Step 1: QC dataset ################ 
sed -e 's/South Asia/SouthAsia/g' -e 's/East Asia/EastAsia/g' -e 's/Middle-East/MiddleEast/g' /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/info.txt > $POWER_DIR/Results/myinfo.txt

Rscript $POWER_DIR/clean_ds.r

#QC data
if [ ! -d "$POWER_DIR/Results/allSNPs" ]; then
    mkdir -p "$POWER_DIR/Results/allSNPs"
fi
if [ ! -d "$POWER_DIR/Results/allSNPs/pca_results" ]; then
    mkdir -p "$POWER_DIR/Results/allSNPs/pca_results"
fi
for ANC in Africa EastAsia European MiddleEast SouthAsia American; do
    grep $ANC $POWER_DIR/Results/myinfo.QC0.txt | cut -f1-2 > $POWER_DIR/Results/idtokeep.$ANC.txt
    $PLINK2 --bfile $TGP_PLINK --maf 0.05 --keep $POWER_DIR/Results/idtokeep.$ANC.txt --make-bed --out $POWER_DIR/Results/allSNPs/TGP_PLINK
    $PLINK2 --bfile $POWER_DIR/Results/allSNPs/TGP_PLINK --indep-pairwise 50 10 0.1 --out $POWER_DIR/Results/allSNPs/TGP_PLINK.allSNPs_pruned
    $PLINK2 --bfile $POWER_DIR/Results/allSNPs/TGP_PLINK --extract $POWER_DIR/Results/allSNPs/TGP_PLINK.allSNPs_pruned.prune.in --make-bed --out $POWER_DIR/Results/allSNPs/TGP_PLINK.allSNPs_pruned
    $PLINK2 --bfile $POWER_DIR/Results/allSNPs/TGP_PLINK.allSNPs_pruned --pca --out $POWER_DIR/Results/allSNPs/pca_results/allSNPs0.$ANC
    
    rm $POWER_DIR/Results/allSNPs/TGP_PLINK.*
done

Rscript $POWER_DIR/split_ds.r

for ANC in Africa EastAsia European MiddleEast SouthAsia American; do
    grep $ANC $POWER_DIR/Results/myinfo.split.txt | cut -f3 | sort | uniq | grep -v POP > $POWER_DIR/Results/listpop.$ANC.txt
done

################ Step 2: Run PCA and ADMIXTURE on each fold ################
list=("allSNPs" "GompertsAirwatCfCells" "HnsccImmuneLandscape" "humanPreimplantationEmbryos" "StemCellsInChronicMyeloidLeukemia")
for STUDY in "${list[@]}"; do
    bash $POWER_DIR/power_region.sh $STUDY
done


