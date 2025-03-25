#!/bin/bash

d=/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/0_Download_data
STUDY_DIR="/scratch1/yaojiani/eQTLAutoimmune/pool1"

while read line
do
    SRAID=$line
    sbatch getFASTQ.sh $SRAID
done < $d/SRR_Acc_List.txt