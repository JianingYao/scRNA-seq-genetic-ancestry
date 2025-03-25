#!/bin/bash

set -e

STUDY=$1

TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
INFO="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"
PLINK2="/project/gazal_569/soft/plink1.9/plink"
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region"
STUDY_DIR=$POWER_DIR/Results/${STUDY}
RESULT_DIR=$POWER_DIR/Results/${STUDY}/pca_results
SNP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results/${STUDY}"
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
    if [ "$STUDY" == "allSNPs" ]; then
        echo "allSNPs"
        $PLINK2 --bfile $TGP_PLINK --maf 0.05 --keep $POWER_DIR/Results/idtokeep.$ANC.txt --make-bed --out TGP_PLINK.$ANC.$STUDY
    else
        echo $STUDY
        $PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $POWER_DIR/Results/idtokeep.$ANC.txt --make-bed --out TGP_PLINK.$ANC.$STUDY 
    fi
    $PLINK2 --bfile TGP_PLINK.$ANC.$STUDY --indep-pairwise 50 10 0.1 --out TGP_PLINK.$ANC.${STUDY}_pruned
    $PLINK2 --bfile TGP_PLINK.$ANC.$STUDY --extract TGP_PLINK.$ANC.${STUDY}_pruned.prune.in --make-bed --out TGP_PLINK.$ANC.${STUDY}_pruned
    #PCA ALL
    $PLINK2 --bfile TGP_PLINK.$ANC.${STUDY}_pruned --pca --out $RESULT_DIR/${STUDY}.$ANC
    # #PCA SPLIT
    num_clusters=10
    while read LINE; do
        for i in $(seq 1 $num_clusters)
        do
            awk -v anc=$ANC '(NR>1) {if($5==anc) {print $0} }' $POWER_DIR/Results/myinfo.split.txt | awk -v i=$i -v pop=$LINE '{if ($9==i && $3==pop) {print $1, $2, "."} else {print $1, $2, "REF"}}' > TGP_PLINK.$ANC.$LINE.split.cluster$i
            $PLINK2 --bfile TGP_PLINK.$ANC.${STUDY}_pruned  --pca --within TGP_PLINK.$ANC.$LINE.split.cluster$i --pca-cluster-names REF --out $RESULT_DIR/${STUDY}.$ANC.$LINE.cluster$i
        done
        #  Run PCA clustering - check the number of PCs to keep
        Rscript $POWER_DIR/pca_results.r $ANC $STUDY $LINE
        #  Run admixture
        if [[ $ANC == "Africa" ]];     then NPOP=4; fi
        if [[ $ANC == "American" ]];   then NPOP=4; fi
        if [[ $ANC == "EastAsia" ]];   then NPOP=4; fi
        if [[ $ANC == "European" ]];   then NPOP=6; fi
        if [[ $ANC == "MiddleEast" ]]; then NPOP=4; fi
        if [[ $ANC == "SouthAsia" ]];  then NPOP=4; fi
        for i in $(seq 1 $num_clusters)
        do
            sbatch --mem=50000 -t 0-03:00 --partition=main --wrap="bash $POWER_DIR/admixture.sh $i $ANC $STUDY $NPOP $LINE"
        done
    done < $POWER_DIR/Results/listpop.$ANC.txt
done

rm $RESULT_DIR/allSNPs0.* $RESULT_DIR/${STUDY}.*.*.cluster* $RESULT_DIR/${STUDY}.* $STUDY_DIR/TGP_PLINK.*.split.cluster*