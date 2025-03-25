#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for DONOR_DIR in /scratch1/yaojiani/processed/eQTLAutoimmune/02_AddRG/*
do 
    DONOR="$(basename $DONOR_DIR)"
    ls $DONOR_DIR/$DONOR*.step2.bam > $DONOR_DIR/$DONOR.txt
    FILE_LIST=$DONOR_DIR/$DONOR.txt 
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/03_mergeBAM/03_merge_eQTLAutoimmune.sh $DONOR $FILE_LIST
done