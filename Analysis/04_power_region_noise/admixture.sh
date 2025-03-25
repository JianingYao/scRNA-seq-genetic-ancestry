#!/bin/bash

set -e

i=$1
ANC=$2
STUDY=$3
NPOP=$4
REGION=$5
POWER_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region_noise"
PLINK2="/project/gazal_569/soft/plink1.9/plink"
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.3.2


awk -v anc=$ANC '(NR>1) {if($5==anc) {print $0} }' $POWER_DIR/myinfo.split.txt | awk -v i=$i -v pop=${REGION} '{if ($9==i && $3==pop) {print "."} else {print $3}}' > $ANC.${REGION}.${STUDY}.cluster${i}.tmp.pop
/project/gazal_569/jianing/sc-RNA_ancestry/bin/admixture_linux-1.3.0/admixture --supervised $ANC.${REGION}.${STUDY}.cluster${i}.tmp.bed $NPOP
paste $ANC.${REGION}.${STUDY}.cluster${i}.tmp.pop $ANC.${REGION}.${STUDY}.cluster${i}.tmp.$NPOP.Q | awk '{if($1==".") {print $2, $3, $4, $5, $6, $7, $8} }' > admixture_results/$STUDY.$ANC.${REGION}.$i.txt

rm $ANC.${REGION}.${STUDY}.cluster${i}.tmp.*
