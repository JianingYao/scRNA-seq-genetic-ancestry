#!/bin/bash

set -e

STUDY=$1
ERROR=$2
ANC=$3

TGP_PLINK="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19"
INFO="/project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc.info"
PLINK2="/project/gazal_569/soft/plink1.9/plink"
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.3.2

POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise"
STUDY_DIR=$POWER_DIR/Results_$ERROR/${STUDY}
RESULT_DIR=$POWER_DIR/Results_$ERROR/${STUDY}/pca_results
SNP_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results/${STUDY}"
cd $STUDY_DIR 

echo $STUDY
echo $ANC
$PLINK2 --bfile $TGP_PLINK --extract $SNP_DIR/$STUDY.TGP.bim --keep $POWER_DIR/idtokeep.$ANC.txt --make-bed --out TGP_PLINK.$ANC.$STUDY
$PLINK2 --bfile TGP_PLINK.$ANC.$STUDY --indep-pairwise 50 10 0.1 --out TGP_PLINK.$ANC.${STUDY}_pruned
$PLINK2 --bfile TGP_PLINK.$ANC.$STUDY --extract TGP_PLINK.$ANC.${STUDY}_pruned.prune.in --make-bed --out TGP_PLINK.$ANC.${STUDY}_pruned
#PCA ALL
# $PLINK2 --bfile TGP_PLINK.$ANC.${STUDY}_pruned --pca --out $RESULT_DIR/${STUDY}.$ANC
#PCA SPLIT
num_clusters=10
while read LINE; do
    for i in $(seq 1 $num_clusters)
    do
        awk -v anc=$ANC '(NR>1) {if($5==anc) {print $0} }' $POWER_DIR/myinfo.split.txt | awk -v i=$i -v pop=$LINE '{if ($9==i && $3==pop) {print $1, $2, "."} else {print $1, $2, "REF"}}' > TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster
        ##a) add noise
        grep REF TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster > TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster.list
        $PLINK2 --bfile TGP_PLINK.$ANC.${STUDY}_pruned --keep   TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster.list --make-bed --out $ANC.$LINE.${STUDY}.tmp1
        $PLINK2 --bfile TGP_PLINK.$ANC.${STUDY}_pruned --remove TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster.list --recode transpose --out $ANC.$LINE.${STUDY}.tmp2
        Rscript $POWER_DIR/add_noise.r $ANC.$LINE.${STUDY} $ERROR
        $PLINK2 --tfile $ANC.$LINE.${STUDY}.tmp2 --make-bed --out $ANC.$LINE.${STUDY}.tmp2
        $PLINK2 --bfile $ANC.$LINE.${STUDY}.tmp1 --bmerge $ANC.$LINE.${STUDY}.tmp2 --make-bed --out $ANC.$LINE.${STUDY}.cluster${i}.tmp
        rm $ANC.$LINE.${STUDY}.tmp?.* TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster.list
        #
        $PLINK2 --bfile $ANC.$LINE.${STUDY}.cluster${i}.tmp  --pca --within TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster --pca-cluster-names REF --out $RESULT_DIR/${STUDY}.$ANC.$LINE.cluster$i
    done
    #Run PCA clustering - check the number of PCs to keep
    Rscript $POWER_DIR/pca_results.r $ANC $STUDY $LINE
    #Run admixture
    if [[ $ANC == "Africa" ]];     then NPOP=4; fi
    if [[ $ANC == "American" ]];   then NPOP=4; fi
    if [[ $ANC == "EastAsia" ]];   then NPOP=4; fi
    if [[ $ANC == "European" ]];   then NPOP=6; fi
    if [[ $ANC == "MiddleEast" ]]; then NPOP=4; fi
    if [[ $ANC == "SouthAsia" ]];  then NPOP=4; fi
    for i in $(seq 1 $num_clusters)
    do
        bash $POWER_DIR/admixture.sh $i $ANC $STUDY $NPOP $LINE
    done
    rm $STUDY_DIR/TGP_PLINK.$ANC.$LINE.${STUDY}.split.cluster $RESULT_DIR/${STUDY}.$ANC.$LINE.*
done < $POWER_DIR/listpop.$ANC.txt


