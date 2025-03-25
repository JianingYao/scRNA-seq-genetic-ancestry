#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for DONOR_DIR in /scratch1/yaojiani/processed/eQTLAutoimmune/02_AddRG/*
do 
    DONOR="$(basename $DONOR_DIR)"
    for FILE in /scratch1/yaojiani/processed/eQTLAutoimmune/02_AddRG/$DONOR/*
    do
        SRA="$(basename $FILE | awk -F'.bam' '{print $1}' | awk -F'_' '{print $3}')"
        sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/02_AddRG/part2_addRG/02_AddRG_eQTLAutoimmune.sh $FILE $DONOR $SRA
    done
done