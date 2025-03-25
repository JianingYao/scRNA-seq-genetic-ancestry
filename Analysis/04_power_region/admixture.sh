#!/bin/bash

set -e

i=$1
ANC=$2
SNP=$3
NPOP=$4
REGION=$5
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region"

awk -v anc=$ANC '(NR>1) {if($5==anc) {print $0} }' $POWER_DIR/Results/myinfo.split.txt | awk -v i=$i -v pop=$REGION '{if ($9==i && $3==pop) {print "."} else {print $3}}' > TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.pop
cp TGP_PLINK.$ANC.${SNP}_pruned.bed  TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.bed
cp TGP_PLINK.$ANC.${SNP}_pruned.bim  TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.bim
cp TGP_PLINK.$ANC.${SNP}_pruned.fam  TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.fam
/project/gazal_569/jianing/sc-RNA_ancestry/bin/admixture_linux-1.3.0/admixture --supervised TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.bed $NPOP
paste TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.pop TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.$NPOP.Q | awk '{if($1==".") {print $2, $3, $4, $5, $6, $7, $8} }' > admixture_results/$SNP.$ANC.${REGION}.$i.txt

rm TGP_PLINK.$ANC.${REGION}.${SNP}.${i}_pruned.*
