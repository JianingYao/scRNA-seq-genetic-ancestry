#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch1/yaojiani/processed/HnsccImmuneLandscape/05_SplitNCigarReads/complete/*.step5.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step5.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/07_ApplyBQSR/07_ApplyBQSR_HnsccImmuneLandscape.sh $FILE_NAME
done